// Copyright (C) 2018(-2024) Reinhard <octave-user@a1.net>
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; If not, see <http://www.gnu.org/licenses/>.

#include "config.h"

#include <algorithm>
#include <atomic>
#include <array>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <memory>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

using namespace std::string_literals;

#if HAVE_NLOPT == 1
#include <nlopt.h>
#endif

#ifdef DEBUG
#define OCTAVE_ENABLE_BOUNDS_CHECK
#define FEM_ASSERT(expr) if (!(expr)) {                 \
          std::cerr << "assertion " << #expr            \
                    << " failed at " << __FILE__ << ":" \
                    << __LINE__ << ":"                  \
                    << __FUNCTION__ << std::endl;       \
          asm("int3");                                  \
          abort();                                      \
     }
#else
#define FEM_ASSERT(expr) static_cast<void>(0)
#endif

#if DEBUG > 1
#define FEM_TRACE(expr) static_cast<void>(std::cout << expr)
#else
#define FEM_TRACE(expr) static_cast<void>(0)
#endif

#include <octave/oct.h>
#include <octave/oct-map.h>
#include <octave/parse.h>

#ifdef DEBUG
#define xelem checkelem
#endif

template <typename T, typename A = std::allocator<T> >
class vector: public std::vector<T, A> {
     typedef std::vector<T, A> base_vector;

public:
     typedef typename base_vector::size_type size_type;
     typedef typename base_vector::allocator_type allocator_type;

     vector() {
     }

     explicit vector(size_type n, const allocator_type& a = allocator_type())
          :base_vector(n, a) {
     }

     T& operator[](size_type i) {
          FEM_ASSERT(i < this->size());
          return base_vector::operator[](i);
     }

     const T& operator[](size_type i) const {
          FEM_ASSERT(i < this->size());
          return base_vector::operator[](i);
     }
};

template <typename T, std::size_t N>
struct array
{
     typedef typename std::array<T, N>::reference reference;
     typedef typename std::array<T, N>::const_reference const_reference;
     typedef typename std::array<T, N>::size_type size_type;
     typedef typename std::array<T, N>::iterator iterator;
     typedef typename std::array<T, N>::const_iterator const_iterator;

     reference operator[](size_type i) {
          FEM_ASSERT(i < data.size());
          return data[i];
     }

     const_reference operator[](size_type i) const {
          FEM_ASSERT(i < data.size());
          return data[i];
     }

     size_type size() const {
          return data.size();
     }

     iterator begin() {
          return data.begin();
     }

     iterator end() {
          return data.end();
     }

     const_iterator begin() const {
          return data.begin();
     }

     const_iterator end() const {
          return data.end();
     }

     std::array<T, N> data;
};

template <typename T>
inline void atomic_fetch_add_float(volatile T& x, const T& dx) noexcept {
     T xprev, xnew;

     __atomic_load(&x, &xprev, __ATOMIC_RELAXED);

     do {
          xnew = xprev + dx;
     } while (!__atomic_compare_exchange(&x, &xprev, &xnew, true, __ATOMIC_SEQ_CST, __ATOMIC_RELAXED));
}

template <typename T>
inline void atomic_fetch_add_float(volatile std::complex<T>& x, const std::complex<T>& dx) noexcept {
     constexpr octave_idx_type N = 2;

     static_assert(sizeof(std::complex<T>) == N * sizeof(T));

     // FIXME: x.real() and x.imag() do not return L-values
     volatile T* const px = reinterpret_cast<volatile double*>(&x);
     volatile T* const pdx = reinterpret_cast<const double*>(&dx);

     for (octave_idx_type i = 0; i < N; ++i) {
          atomic_fetch_add_float(px[i], pdx[i]);
     }
}

template <typename T>
struct PostProcTypeTraits;

template <>
struct PostProcTypeTraits<double> {
     typedef Matrix MatrixType;
     typedef ColumnVector ColumnVectorType;
     typedef RowVector RowVectorType;
     typedef NDArray NDArrayType;

     static double
     EffectiveAmplitude(double x, double y) {
          // return the instantaneous value of x * y
          return x * y;
     }
};

template <>
struct PostProcTypeTraits<std::complex<double> > {
     typedef ComplexMatrix MatrixType;
     typedef ComplexColumnVector ColumnVectorType;
     typedef ComplexRowVector RowVectorType;
     typedef ComplexNDArray NDArrayType;

     static double
     EffectiveAmplitude(const std::complex<double>& x, const std::complex<double>& y) {
          // return the effective amplitude of x * y
          // where x and y are complex amplitudes of harmonic quantities
          // e.g. x(t) = realpart(x * exp(j * omega * t))
          //      y(t) = realpart(y * exp(j * omega * t))
          return std::abs(x) * std::abs(y) / 2.;
     }
};

class DofMap {
public:
     enum DomainType: unsigned {
          DO_STRUCTURAL   = 0x01000000u,
          DO_THERMAL      = 0x02000000u,
          DO_ACOUSTICS    = 0x04000000u,
          DO_FLUID_STRUCT = 0x08000000u
     };

     enum: unsigned {
          DO_MASK       = 0xFF000000u
     };

     enum NodalDofType: octave_idx_type {
          NDOF_DISPLACEMENT,
          NDOF_TEMPERATURE,
          NDOF_VELOCITY_POT,
          NDOF_COUNT
     };

     enum ElementType {
          ELEM_RBE3 = 0,
          ELEM_JOINT,
          ELEM_TYPE_COUNT,
          ELEM_NODOF = -1
     };

     explicit DofMap(DomainType eDomain, const int32NDArray& ndof, const array<int32NDArray, ELEM_TYPE_COUNT>& edof, octave_idx_type totdof)
          :eDomain(eDomain), ndof(ndof), edof(edof), totdof(totdof) {

          std::fill(std::begin(ioffset), std::end(ioffset), INVALID_OFFSET);

          switch (eDomain) {
          case DO_STRUCTURAL:
               ioffset[NDOF_DISPLACEMENT] = 0L;
               break;

          case DO_THERMAL:
               ioffset[NDOF_TEMPERATURE] = 0L;
               break;

          case DO_ACOUSTICS:
               ioffset[NDOF_VELOCITY_POT] = 0L;
               break;

          case DO_FLUID_STRUCT:
               ioffset[NDOF_DISPLACEMENT] = 0L;
               ioffset[NDOF_VELOCITY_POT] = iGetNodeMaxDofIndex(DO_STRUCTURAL);
               break;

          default:
               throw std::runtime_error("dof_map: invalid value for dof_map.domain");
          }

#ifdef DEBUG
          for (octave_idx_type i = 0; i < NDOF_COUNT; ++i) {
               FEM_ASSERT(ioffset[i] == INVALID_OFFSET || (ioffset[i] >= 0 && ioffset[i] < ndof.columns()));
          }
#endif
     }

     octave_idx_type GetNodeDofIndex(octave_idx_type inode, NodalDofType etype, octave_idx_type idof) const {
          FEM_ASSERT(ioffset[etype] != INVALID_OFFSET);

          return ndof.xelem(inode + ndof.rows() * (ioffset[etype] + idof));
     }

     octave_idx_type GetElemDofIndex(ElementType eElemType, octave_idx_type ielem, octave_idx_type idof) const {
          const auto& edoftype = edof[eElemType];
          return edoftype.xelem(ielem + edoftype.rows() * idof).value();
     }

     octave_idx_type iGetNumDof() const {
          return totdof;
     }

     octave_idx_type iGetNumNodes() const {
          return ndof.rows();
     }

     octave_idx_type iGetNodeMaxDofIndex() const {
          return ndof.columns();
     }

     static octave_idx_type iGetNodeMaxDofIndex(DomainType eDomain) {
          switch (eDomain) {
          case DO_STRUCTURAL:
               return 6;

          case DO_THERMAL:
          case DO_ACOUSTICS:
               return 1;

          case DO_FLUID_STRUCT:
               return 7;

          default:
               return -1;
          }
     }

     DomainType GetDomain() const {
          return eDomain;
     }

private:
     static constexpr octave_idx_type INVALID_OFFSET = std::numeric_limits<octave_idx_type>::min();

     const DomainType eDomain;
     const int32NDArray ndof;
     array<int32NDArray, ELEM_TYPE_COUNT> edof;
     array<octave_idx_type, NDOF_COUNT> ioffset;
     const octave_idx_type totdof;
};

class ParallelOptions {
public:
     explicit ParallelOptions(const octave_scalar_map& oDof)
          :iNumThreadsAss{1},
           iMultiThreadThreshold{10000} {

          const auto iter_parallel = oDof.seek("parallel");

          if (iter_parallel != oDof.end()) {
               const octave_scalar_map ma_parallel = oDof.contents(iter_parallel).scalar_map_value();

               const auto iter_threads_ass = ma_parallel.seek("threads_ass");

               if (iter_threads_ass != ma_parallel.end()) {
                    iNumThreadsAss = ma_parallel.contents(iter_threads_ass).int_value();
               }

               const auto iter_threshold_elem = ma_parallel.seek("threshold_elem");

               if (iter_threshold_elem != ma_parallel.end()) {
                    iMultiThreadThreshold = ma_parallel.contents(iter_threshold_elem).int_value();
               }
          }
     }

     ParallelOptions(octave_idx_type iNumThreadsAss, octave_idx_type iMultiThreadThreshold)
          :iNumThreadsAss(iNumThreadsAss),
           iMultiThreadThreshold(iMultiThreadThreshold) {
     }

     octave_idx_type iGetNumThreadsAss() const {
          return iNumThreadsAss;
     }

     octave_idx_type iGetMultiThreadThreshold() const {
          return iMultiThreadThreshold;
     }

private:
     octave_idx_type iNumThreadsAss;
     octave_idx_type iMultiThreadThreshold;
};

constexpr octave_idx_type DofMap::INVALID_OFFSET;

class MatrixAss;

class MeshInfo {
public:
     enum InfoType {
          JACOBIAN_DET = 0,
          JACOBIAN_DET_A = 1,
          INFO_COUNT
     };

     enum InfoStat {
          STAT_MIN = 0,
          STAT_MAX,
          STAT_MEAN,
          STAT_COUNT
     };

     MeshInfo() {
          Reset();
     }

     octave_value Get() const {
          octave_scalar_map oMeshInfo;

          static constexpr char szStatName[STAT_COUNT][5] = {
               "min",
               "max",
               "mean"
          };

          static constexpr char szInfoName[INFO_COUNT][6] = {
               "detJ",
               "detJA"
          };

          for (octave_idx_type i = 0; i < INFO_COUNT; ++i) {
               octave_scalar_map oMeshStat;
               const InfoType eInfoType = static_cast<InfoType>(i);

               for (octave_idx_type j = 0; j < STAT_COUNT; ++j) {
                    const InfoStat eInfoStat = static_cast<InfoStat>(j);
                    oMeshStat.assign(szStatName[j], dGet(eInfoType, eInfoStat));
               }

               oMeshStat.assign("count", rgData[i].iCount);
               oMeshInfo.assign(szInfoName[i], oMeshStat);
          }

          return oMeshInfo;
     }

     void Reset() {
          for (auto i = std::begin(rgData); i != std::end(rgData); ++i) {
               i->rgValue[STAT_MIN] = std::numeric_limits<double>::max();
               i->rgValue[STAT_MAX] = -i->rgValue[STAT_MIN];
               i->rgValue[STAT_MEAN] = 0;
               i->iCount = 0;
          }
     }

     void Add(InfoType eInfoType, double dValue) {
          Data& oData = rgData[eInfoType];
          oData.rgValue[STAT_MIN] = std::min(oData.rgValue[STAT_MIN], dValue);
          oData.rgValue[STAT_MAX] = std::max(oData.rgValue[STAT_MAX], dValue);
          oData.rgValue[STAT_MEAN] += dValue;
          ++oData.iCount;
     }

     void Add(const MeshInfo& oMeshInfo) {
          for (octave_idx_type i = 0; i < INFO_COUNT; ++i) {
               Data& oData = rgData[i];
               const Data& oDataAdd = oMeshInfo.rgData[i];
               oData.rgValue[STAT_MIN] = std::min(oData.rgValue[STAT_MIN], oDataAdd.rgValue[STAT_MIN]);
               oData.rgValue[STAT_MAX] = std::max(oData.rgValue[STAT_MAX], oDataAdd.rgValue[STAT_MAX]);
               oData.rgValue[STAT_MEAN] += oDataAdd.rgValue[STAT_MEAN];
               oData.iCount += oDataAdd.iCount;
          }
     }

     double dGet(InfoType eInfoType, InfoStat eInfoStat) const {
          double dValue = rgData[eInfoType].rgValue[eInfoStat];

          if (eInfoStat == STAT_MEAN) {
               dValue /= rgData[eInfoType].iCount;
          }

          return dValue;
     }

private:
     struct Data {
          array<double, STAT_COUNT> rgValue;
          int iCount;
     };

     array<Data, INFO_COUNT> rgData;
};

class Material
{
public:
     enum MatType {
          MAT_TYPE_SOLID,
          MAT_TYPE_THERMAL,
          MAT_TYPE_FLUID
     };

     Material(MatType eMatType, const Matrix& C, double rho, double alpha, double beta, double gamma, const Matrix& k, double cp, double c, double eta, double zeta, double tan_delta)
          :eMatType(eMatType),
           rho(rho),
           alpha(alpha),
           beta(beta),
           gamma(gamma),
           cp(cp),
           c(c),
           eta(eta),
           zeta(zeta),
           tan_delta(tan_delta),
           C(C),
           k(k) {

          if (eMatType == MAT_TYPE_SOLID) {
               FEM_ASSERT(C.rows() == 6);
               FEM_ASSERT(C.columns() == 6);

               const double a = C.xelem(0);
               const double b = C.xelem(5 + 6 * 5) / a;

               E = ((4 * b * b - 3 * b) * a) / (b - 1);
               nu = (2 * b - 1) / (2 * b - 2);
          } else {
               E = nu = -1.;
          }
     }

     static vector<Material> ExtractMaterialData(const octave_map& material_data, const enum DofMap::DomainType eDomain) {
          const auto iterE = material_data.seek("E");
          const auto iternu = material_data.seek("nu");
          const auto iterC = material_data.seek("C");
          const auto iterRho = material_data.seek("rho");
          const auto iterAlpha = material_data.seek("alpha");
          const auto iterBeta = material_data.seek("beta");
          const auto iterGamma = material_data.seek("gamma");
          const auto iterk = material_data.seek("k");
          const auto itercp = material_data.seek("cp");
          const auto iterc = material_data.seek("c");
          const auto itereta = material_data.seek("eta");
          const auto iterzeta = material_data.seek("zeta");
          const auto itertan_delta = material_data.seek("tan_delta");

          const Cell emptyCell(material_data.dims());
          const Cell cellC = iterC != material_data.end() ? material_data.contents(iterC) : emptyCell;
          const Cell cellE = iterE != material_data.end() ? material_data.contents(iterE) : emptyCell;
          const Cell cellnu = iternu != material_data.end() ? material_data.contents(iternu) : emptyCell;
          const Cell cellRho = iterRho != material_data.end() ? material_data.contents(iterRho) : emptyCell;
          const Cell cellAlpha = iterAlpha != material_data.end() ? material_data.contents(iterAlpha) : emptyCell;
          const Cell cellBeta = iterBeta != material_data.end() ? material_data.contents(iterBeta) : emptyCell;
          const Cell cellGamma = iterGamma != material_data.end() ? material_data.contents(iterGamma) : emptyCell;
          const Cell cellk = iterk != material_data.end() ? material_data.contents(iterk) : emptyCell;
          const Cell cellcp = itercp != material_data.end() ? material_data.contents(itercp) : emptyCell;
          const Cell cellc = iterc != material_data.end() ? material_data.contents(iterc) : emptyCell;
          const Cell celleta = itereta != material_data.end() ? material_data.contents(itereta) : emptyCell;
          const Cell cellzeta = iterzeta != material_data.end() ? material_data.contents(iterzeta) : emptyCell;
          const Cell celltan_delta = itertan_delta != material_data.end() ? material_data.contents(itertan_delta) : emptyCell;

          vector<Material> rgMaterials;

          rgMaterials.reserve(material_data.numel());

          Matrix C;

          for (octave_idx_type i = 0; i < material_data.numel(); ++i) {
               const bool bHaveC = !cellC.xelem(i).isempty();
               const bool bHaveE = !cellE.xelem(i).isempty();
               const bool bHavenu = !cellnu.xelem(i).isempty();
               const bool bHaveRho = !cellRho.xelem(i).isempty();
               const bool bHaveAlpha = !cellAlpha.xelem(i).isempty();
               const bool bHaveBeta = !cellBeta.xelem(i).isempty();
               const bool bHaveGamma = !cellGamma.xelem(i).isempty();
               const bool bHavek = !cellk.xelem(i).isempty();
               const bool bHavecp = !cellcp.xelem(i).isempty();
               const bool bHavec = !cellc.xelem(i).isempty();
               const bool bHaveeta = !celleta.xelem(i).isempty();
               const bool bHavezeta = !cellzeta.xelem(i).isempty();
               const bool bHaveElasticity = bHaveC || (bHaveE && bHavenu);
               const bool bHaveTanDelta = !celltan_delta.xelem(i).isempty();

               Material::MatType eMatType;

               switch (eDomain) {
               case DofMap::DO_STRUCTURAL:
                    eMatType = Material::MAT_TYPE_SOLID;
                    break;

               case DofMap::DO_THERMAL:
                    eMatType = Material::MAT_TYPE_THERMAL;
                    break;

               case DofMap::DO_ACOUSTICS:
                    eMatType = Material::MAT_TYPE_FLUID;
                    break;

               case DofMap::DO_FLUID_STRUCT:
                    if (!bHavec && bHaveElasticity) {
                         eMatType = Material::MAT_TYPE_SOLID;
                    } else if (bHavec && !bHaveElasticity) {
                         eMatType = Material::MAT_TYPE_FLUID;
                    } else {
                         throw std::runtime_error("material: mesh.material_data is not consistent for fluid structure interaction");
                    }
                    break;

               default:
                    throw std::logic_error("material: unsupported value for dof_map.domain");
               }

               switch (eMatType) {
               case Material::MAT_TYPE_SOLID:
               case Material::MAT_TYPE_THERMAL:
                    if (bHavec || bHaveeta || bHavezeta) {
                         throw std::runtime_error("material: fields \"c\", \"eta\" and \"zeta\" are not valid properites for solids in mesh.material_data");
                    }
                    break;
               case Material::MAT_TYPE_FLUID:
                    if (bHaveElasticity || bHavek || bHavecp) {
                         throw std::runtime_error("material: fields \"E\", \"nu\", \"C\", \"k\" and \"cp\" ar not valid properties for fluids in mesh.material_data");
                    }
                    break;
               }

               switch (eMatType) {
               case Material::MAT_TYPE_SOLID:
                    if (bHaveC) {
                         if (bHaveE || bHavenu) {
                              throw std::runtime_error("material: redundant material properties in field mesh.material_data");
                         }

                         C = cellC.xelem(i).matrix_value();

#if OCTAVE_MAJOR_VERSION < 6
                         if (error_state) {
                              throw std::runtime_error("material: mesh.material_data.C must be matrix");
                         }
#endif
                    } else {
                         if (!(bHaveE && bHavenu)) {
                              throw std::runtime_error("material: fields \"E\" and \"nu\" not found in mesh.material_data");
                         }

                         const double E = cellE.xelem(i).scalar_value();

#if OCTAVE_MAJOR_VERSION < 6
                         if (error_state) {
                              throw std::runtime_error("material: mesh.material_data.E must be a real scalar");
                         }
#endif

                         if (E <= 0.) {
                              throw std::runtime_error("material: mesh.material_data.E must be greater than zero");
                         }

                         const double nu = cellnu.xelem(i).scalar_value();

#if OCTAVE_MAJOR_VERSION < 6
                         if (error_state) {
                              throw std::runtime_error("material: \" in mesh.material_data must be a real scalar");
                         }
#endif
                         if (nu < 0. || nu >= 0.5) {
                              throw std::runtime_error("material: field \"nu\" in mesh.material_data must be between 0 and 0.5");
                         }

                         C = Material::IsotropicElasticity(E, nu);
                    }

                    if (C.rows() != 6 || C.columns() != 6) {
                         throw std::runtime_error("material: size of constitutive matrix mesh.material_data.C is not valid in argument mesh");
                    }

                    if (!C.issymmetric()) {
                         throw std::runtime_error("material: mesh.material_data.C is not symmetric");
                    }

                    for (octave_idx_type j = 0; j < C.rows(); ++j) {
                         for (octave_idx_type k = 0; k < C.columns(); ++k) {
                              if (!std::isfinite(C.xelem(j, k))) {
                                   throw std::runtime_error("material: mesh.material_data.C is not finite");
                              }
                         }
                    }

                    if (!bHaveRho) {
                         throw std::runtime_error("material: missing field \"rho\" in mesh.material_data");
                    }
                    break;

               case Material::MAT_TYPE_THERMAL:
                    if (!(bHavek && bHavecp && bHaveRho)) {
                         throw std::runtime_error("material: missing fields mesh.material_data.k, mesh.material_data.cp and mesh.material_data.rho");
                    }
                    break;

               case Material::MAT_TYPE_FLUID:
                    if (!(bHavec && bHaveRho)) {
                         throw std::runtime_error("material: missing field mesh.material_data.c and mesh.material_data.rho");
                    }
                    break;

               default:
                    throw std::logic_error("material: unsupported material type");
               }


               const double rho = cellRho.xelem(i).scalar_value();

#if OCTAVE_MAJOR_VERSION < 6
               if (error_state) {
                    throw std::runtime_error("material: mesh.material_data.rho is not a valid scalar in argument mesh");
               }
#endif
               if (rho < 0.) {
                    throw std::runtime_error("material: mesh.material_data.rho must be greater than or equal to zero");
               }

               const double alpha = bHaveAlpha ? cellAlpha.xelem(i).scalar_value() : 0.;

#if OCTAVE_MAJOR_VERSION < 6
               if (error_state) {
                    throw std::runtime_error("material: mesh.material_data.alpha is not a valid scalar in argument mesh");
               }
#endif
               if (bHaveAlpha && alpha < 0.) {
                    throw std::runtime_error("material: mesh.material_data.alpha must be greater than or equal to zero");
               }

               const double beta = bHaveBeta ? cellBeta.xelem(i).scalar_value() : 0.;

#if OCTAVE_MAJOR_VERSION < 6
               if (error_state) {
                    throw std::runtime_error("material: mesh.material_data.beta is not a valid scalar in argument mesh");
               }
#endif
               if (bHaveBeta && beta < 0.) {
                    throw std::runtime_error("material: mesh.material_data.beta must be greater than or equal to zero");
               }

               const double gamma = bHaveGamma ? cellGamma.xelem(i).scalar_value() : 0.;

               const Matrix k = bHavek ? cellk.xelem(i).matrix_value() : Matrix(3, 3, 0.);

               if (k.rows() != 3 || k.columns() != 3 || !k.issymmetric()) {
                    throw std::runtime_error("material: mesh.material_data.k must be a real symmetric 3x3 matrix");
               }

               for (octave_idx_type l = 0; l < k.rows(); ++l) {
                    for (octave_idx_type m = 0; m < k.columns(); ++m) {
                         if (!std::isfinite(k.xelem(l, m))) {
                              throw std::runtime_error("material: mesh.material_data.k is not finite");
                         }
                    }
               }

               const double cp = bHavecp ? cellcp.xelem(i).scalar_value() : 0.;

               if (bHavecp && cp <= 0.) {
                    throw std::runtime_error("material: mesh.material_data.cp must be greater than zero");
               }

               const double c = bHavec ? cellc.xelem(i).scalar_value() : 0.;

               if (bHavec && c <= 0.) {
                    throw std::runtime_error("material: mesh.material_data.c must be greater than zero");
               }

               const double eta = bHaveeta ? celleta.xelem(i).scalar_value() : 0.;

               if (bHaveeta && eta < 0.) {
                    throw std::runtime_error("material: mesh.material_data.eta must be greater than or equal to zero");
               }

               const double zeta = bHavezeta ? cellzeta.xelem(i).scalar_value() : 0.;

               if (bHavezeta && zeta < 0.) {
                    throw std::runtime_error("material: mesh.material_data.zeta must be greater than or equal to zero");
               }

               const double tan_delta = bHaveTanDelta ? celltan_delta.xelem(i).scalar_value() : 0.;

               if (tan_delta < 0.) {
                    throw std::runtime_error("material: mesh.material_data.tan_delta must be greater than or equal to zero");
               }

               rgMaterials.emplace_back(eMatType, C, rho, alpha, beta, gamma, k, cp, c, eta, zeta, tan_delta);
          }

          return rgMaterials;
     }

     MatType GetMaterialType() const {
          return eMatType;
     }

     const Matrix& LinearElasticity() const {
          return C;
     }

     const Matrix& ThermalConductivity() const {
          return k;
     }

     double HeatCapacity() const {
          return cp;
     }

     double ThermalExpansion() const {
          return gamma;
     }

     double YoungsModulus() const {
          return E;
     }

     double PoissonsRatio() const {
          return nu;
     }

     double ShearModulus() const {
          return E / (2 * (1 + nu));
     }

     double Density() const { return rho; }
     double AlphaDamping() const { return alpha; }
     double BetaDamping() const { return beta; }
     double TanDeltaDamping() const { return tan_delta; }
     double SpeedOfSound() const { return c; }
     double ShearViscosity() const { return eta; }
     double VolumeViscosity() const { return zeta; }

private:
     static Matrix IsotropicElasticity(const double E, const double nu) {
          Matrix C(6, 6, 0.);

          const double a = nu / (1 - nu);
          const double b = (1 - 2 * nu) / (2 * (1 - nu));
          const double d = E * (1 - nu) / ((1 + nu) * (1 - 2 * nu));

          for (octave_idx_type j = 0; j < 3; ++j) {
               for (octave_idx_type i = 0; i < 3; ++i) {
                    C.xelem(i + 6 * j) = (i == j) ? d : a * d;
               }

               C.xelem(j + 3 + 6 * (j + 3)) = b * d;
          }

          return C;
     }

     MatType eMatType;
     double E, nu, rho, alpha, beta, gamma, cp, c, eta, zeta, tan_delta;
     Matrix C, k;
};

class IntegrationRule
{
public:
     IntegrationRule()
          :iNumEvalPoints{0}, iNumDirections{0} {
          FEM_ASSERT(bInvariant());
     }

     void SetNumEvalPoints(octave_idx_type iNumEvalPnt, octave_idx_type iNumDir) {
          FEM_ASSERT(bInvariant());
          r.resize(iNumEvalPnt, iNumDir);
          alpha.resize(iNumEvalPnt);
          iNumEvalPoints = iNumEvalPnt;
          iNumDirections = iNumDir;
          FEM_ASSERT(bInvariant());
     }

     void SetPosition(octave_idx_type iEvalPnt, octave_idx_type iDirection, double ri) {
          FEM_ASSERT(bInvariant());
          r.xelem(iEvalPnt + iNumEvalPoints * iDirection) = ri;
          FEM_ASSERT(bInvariant());
     }

     void SetWeight(octave_idx_type iEvalPnt, double alphai) {
          FEM_ASSERT(bInvariant());
          alpha.xelem(iEvalPnt) = alphai;
          FEM_ASSERT(bInvariant());
     }

     double dGetPosition(octave_idx_type iEvalPnt, octave_idx_type iDirection) const {
          FEM_ASSERT(bInvariant());
          return r.xelem(iEvalPnt + iNumEvalPoints * iDirection);
     }

     double dGetWeight(octave_idx_type iEvalPnt) const {
          FEM_ASSERT(bInvariant());
          return alpha.xelem(iEvalPnt);
     }

     octave_idx_type iGetNumDirections() const {
          FEM_ASSERT(bInvariant());
          return iNumDirections;
     }

     octave_idx_type iGetNumEvalPoints() const {
          FEM_ASSERT(bInvariant());
          return iNumEvalPoints;
     }

private:
#ifdef DEBUG
     bool bInvariant() const {
          FEM_ASSERT(alpha.rows() == iNumEvalPoints);
          FEM_ASSERT(r.rows() == iNumEvalPoints);
          FEM_ASSERT(r.columns() == iNumDirections);

          return true;
     }
#endif
     Matrix r;
     octave_idx_type iNumEvalPoints, iNumDirections;
     ColumnVector alpha;
};

struct StrainField {
     StrainField(const octave_map& load_case, const Matrix& nodes) {
          const auto iterTemperature = load_case.seek("dTheta");

          if (iterTemperature != load_case.end()) {
               rgTemperature = load_case.contents(iterTemperature);
          }

          const auto iterRefStrain = load_case.seek("epsilon0");

          if (iterRefStrain != load_case.end()) {
               rgRefStrain = load_case.contents(iterRefStrain);
          }

          for (octave_idx_type i = 0; i < rgTemperature.numel(); ++i) {
               const octave_value ov_dTheta = rgTemperature.xelem(i);

               if (!(ov_dTheta.is_matrix_type() &&
                     ov_dTheta.OV_ISREAL() &&
                     ov_dTheta.rows() == nodes.rows() &&
                     ov_dTheta.columns() == 1)) {
                    throw std::runtime_error("strain field: argument load_case.dTheta must be a real column vector with the same number of rows like mesh.nodes");
               }
          }

          for (octave_idx_type i = 0; i < rgRefStrain.numel(); ++i) {
               const octave_value ov_RefStrain = rgRefStrain.xelem(i);

               if (!(ov_RefStrain.isstruct() && ov_RefStrain.numel() == 1)) {
                    throw std::runtime_error("strain field: argument load_case.epsilon0 must be a scalar struct");
               }
          }
     }

     Cell rgTemperature;
     Cell rgRefStrain;
};

struct PreStressField {
     PreStressField(const octave_map& load_case, const Matrix& nodes) {
          const auto iterRefStress = load_case.seek("tau0");

          if (iterRefStress != load_case.end()) {
               rgPreStress = load_case.contents(iterRefStress);
          }

          for (octave_idx_type i = 0; i < rgPreStress.numel(); ++i) {
               const octave_value ov_RefStress = rgPreStress.xelem(i);

               if (!(ov_RefStress.isstruct() && ov_RefStress.numel() == 1)) {
                    throw std::runtime_error("stress field: argument load_case.tau0 must be a scalar struct");
               }
          }
     }

     Cell rgPreStress;
};

struct GravityLoad {
     explicit GravityLoad(const octave_map& load_case)
          :g(3, 0) {
          const auto iterg = load_case.seek("g");

          if (iterg != load_case.end()) {
               const Cell cellg = load_case.contents(iterg);

               FEM_ASSERT(cellg.numel() == load_case.numel());

               constexpr octave_idx_type grows = 3;
               g.resize(grows, load_case.numel(), 0.);

               for (octave_idx_type j = 0; j < load_case.numel(); ++j) {
                    const octave_value ovg = cellg.xelem(j);

                    if (!(ovg.isreal() && ovg.is_matrix_type() && ovg.rows() == grows && ovg.columns() == 1)) {
                         throw std::runtime_error("gravity load: load_case.g must be a real 3x1 vector");
                    }

                    const ColumnVector gj = ovg.column_vector_value();

                    for (octave_idx_type i = 0; i < grows; ++i) {
                         g.xelem(i + grows * j) = gj.xelem(i);
                    }
               }
          }
     }

     Matrix g;
};

struct CentrifugalLoad {
     explicit CentrifugalLoad(const octave_map& load_case)
          :WxWx(dim_vector(3, 3, 0)), Wx(dim_vector(3, 3, 0)) {
          const auto iterW = load_case.seek("omega");
          const auto iterWP = load_case.seek("omegadot");
          const auto iterWq = load_case.seek("omegaq");

          if (iterW != load_case.end() && iterWq != load_case.end()) {
               throw std::runtime_error("centrifugal load: only one parameter out of load_case.omega and load_case.omegaq may be provided");
          }

          if (iterW != load_case.end()) {
               const Cell cellW = load_case.contents(iterW);

               FEM_ASSERT(cellW.numel() == load_case.numel());

               WxWx.resize(dim_vector(3, 3, load_case.numel()), 0.);
               Wx.resize(dim_vector(3, 3, load_case.numel()), 0.);

               for (octave_idx_type j = 0; j < load_case.numel(); ++j) {
                    const octave_value ovW = cellW.xelem(j);

                    if (!(ovW.isreal() && ovW.is_matrix_type() && ovW.rows() == 3 && ovW.columns() == 1)) {
                         throw std::runtime_error("centrifugal load: load_case.omega must be a real 3x1 vector");
                    }

                    const ColumnVector Wj = ovW.column_vector_value();

                    SkewWDotSkewW(WxWx, j, Wj);
                    SkewW(Wx, j, Wj);
               }
          } else if (iterWq != load_case.end()) {
               const Cell cellWq = load_case.contents(iterWq);

               FEM_ASSERT(cellWq.numel() == load_case.numel());

               WxWx.resize(dim_vector(3, 3, load_case.numel()), 0.);

               for (octave_idx_type j = 0; j < load_case.numel(); ++j) {
                    const octave_value ovWq = cellWq.xelem(j);

                    if (!(ovWq.isreal() && ovWq.is_matrix_type() && ovWq.rows() == 6 && ovWq.columns() == 1)) {
                         throw std::runtime_error("centrifugal load: load_case.omegaq must be a real 6x1 vector");
                    }

                    const ColumnVector Wqj = ovWq.column_vector_value();

                    SkewWqSkewWq(WxWx, j, Wqj);
               }
          }

          if (iterWP != load_case.end()) {
               const Cell cellWP = load_case.contents(iterWP);

               FEM_ASSERT(cellWP.numel() == load_case.numel());

               WPx.resize(dim_vector(3, 3, load_case.numel()), 0.);

               for (octave_idx_type j = 0; j < load_case.numel(); ++j) {
                    const octave_value ovWP = cellWP.xelem(j);

                    if (!(ovWP.isreal() && ovWP.is_matrix_type() && ovWP.rows() == 3 && ovWP.columns() == 1)) {
                         throw std::runtime_error("angular acceleration load: load_case.omegaP must be a real 3x1 vector");
                    }

                    const ColumnVector WPj = ovWP.column_vector_value();

                    SkewW(WPx, j, WPj);
               }
          }
     }

     static void SkewWDotSkewW(NDArray& WxWx, octave_idx_type j, const ColumnVector& Wj) {
          // WxWx = \tilde{\omega} \, \tilde{\omega}

          // Wj = \begin{pmatrix}
          //      \omega_1
          //      \omega_2
          //      \omega_3
          //      \end{pmatrix}

          WxWx.xelem(0, 0, j) = -std::pow(Wj.xelem(2), 2) - std::pow(Wj.xelem(1), 2);
          WxWx.xelem(0, 1, j) = Wj.xelem(0) * Wj.xelem(1);
          WxWx.xelem(0, 2, j) = Wj.xelem(0) * Wj.xelem(2);
          WxWx.xelem(1, 0, j) = Wj.xelem(0) * Wj.xelem(1);
          WxWx.xelem(1, 1, j) = -std::pow(Wj.xelem(2), 2) - std::pow(Wj.xelem(0), 2);
          WxWx.xelem(1, 2, j) = Wj.xelem(1) * Wj.xelem(2);
          WxWx.xelem(2, 0, j) = Wj.xelem(0) * Wj.xelem(2);
          WxWx.xelem(2, 1, j) = Wj.xelem(1) * Wj.xelem(2);
          WxWx.xelem(2, 2, j) = -std::pow(Wj.xelem(1), 2) - std::pow(Wj.xelem(0), 2);
     }

     static void SkewWqSkewWq(NDArray& WxWx, octave_idx_type j, const ColumnVector& Wqj) {
          // WxWx = \tilde{\omega} \, \tilde{\omega}

          // Wqj = \begin{pmatrix}
          //       \omega_1^2
          //       \omega_2^2
          //       \omega_3^2
          //       \omega_1 \, omega_2
          //       \omega_2 \, \omega_3
          //       \omega_3 \, \omega_1
          //       \end{pmatrix}

          WxWx.xelem(0, 0, j) = -Wqj.xelem(2) - Wqj.xelem(1);
          WxWx.xelem(0, 1, j) = Wqj.xelem(3);
          WxWx.xelem(0, 2, j) = Wqj.xelem(5);
          WxWx.xelem(1, 0, j) = Wqj.xelem(3);
          WxWx.xelem(1, 1, j) = -Wqj.xelem(2) - Wqj.xelem(0);
          WxWx.xelem(1, 2, j) = Wqj.xelem(4);
          WxWx.xelem(2, 0, j) = Wqj.xelem(5);
          WxWx.xelem(2, 1, j) = Wqj.xelem(4);
          WxWx.xelem(2, 2, j) = -Wqj.xelem(1) - Wqj.xelem(0);
     }

     static void SkewW(NDArray& Wx, octave_idx_type j, const ColumnVector& Wj) {
          Wx.xelem(0, 0, j) = 0;
          Wx.xelem(0, 1, j) = -Wj.xelem(2);
          Wx.xelem(0, 2, j) = Wj.xelem(1);
          Wx.xelem(1, 0, j) = Wj.xelem(2);
          Wx.xelem(1, 1, j) = 0;
          Wx.xelem(1, 2, j) = -Wj.xelem(0);
          Wx.xelem(2, 0, j) = -Wj.xelem(1);
          Wx.xelem(2, 1, j) = Wj.xelem(0);
          Wx.xelem(2, 2, j) = 0;
     }

     NDArray WxWx, Wx, WPx;
};

class PostProcData;

class ElementTypes {
public:
     enum TypeId {
          ELEM_ISO8 = 0,
          ELEM_ISO20,
          ELEM_ISO20R,
          ELEM_ISO27,
          ELEM_PENTA15,
          ELEM_TET10H,
          ELEM_TET10,
          ELEM_TET20,
          ELEM_BEAM2,
          ELEM_RBE3,
          ELEM_JOINT,
          ELEM_SPRING,
          ELEM_DASHPOT,
          ELEM_RIGID_BODY,
          ELEM_SFNCON4,
          ELEM_SFNCON6,
          ELEM_SFNCON6H,
          ELEM_SFNCON8,
          ELEM_SFNCON8R,
          ELEM_SFNCON9,
          ELEM_SFNCON10,
          ELEM_SFNCON4S,
          ELEM_SFNCON6S,
          ELEM_SFNCON6HS,
          ELEM_SFNCON8S,
          ELEM_SFNCON8RS,
          ELEM_SFNCON9S,
          ELEM_SFNCON10S,
          ELEM_PRESSURE_ISO4,
          ELEM_PRESSURE_QUAD8,
          ELEM_PRESSURE_QUAD8R,
          ELEM_PRESSURE_QUAD9,
          ELEM_PRESSURE_TRIA6,
          ELEM_PRESSURE_TRIA6H,
          ELEM_PRESSURE_TRIA10,
          ELEM_STRUCT_FORCE,
          ELEM_THERM_CONV_ISO4,
          ELEM_THERM_CONV_QUAD8,
          ELEM_THERM_CONV_QUAD8R,
          ELEM_THERM_CONV_QUAD9,
          ELEM_THERM_CONV_TRIA6,
          ELEM_THERM_CONV_TRIA6H,
          ELEM_THERM_CONV_TRIA10,
          ELEM_THERM_CONSTR,
          ELEM_HEAT_SOURCE_ISO4,
          ELEM_HEAT_SOURCE_QUAD8,
          ELEM_HEAT_SOURCE_QUAD8R,
          ELEM_HEAT_SOURCE_QUAD9,
          ELEM_HEAT_SOURCE_TRIA6,
          ELEM_HEAT_SOURCE_TRIA6H,
          ELEM_HEAT_SOURCE_TRIA10,
          ELEM_PARTICLE_VEL_ISO4,
          ELEM_PARTICLE_VEL_QUAD8,
          ELEM_PARTICLE_VEL_QUAD8R,
          ELEM_PARTICLE_VEL_QUAD9,
          ELEM_PARTICLE_VEL_TRIA6,
          ELEM_PARTICLE_VEL_TRIA6H,
          ELEM_PARTICLE_VEL_TRIA10,
          ELEM_ACOUSTIC_IMPE_ISO4,
          ELEM_ACOUSTIC_IMPE_QUAD8,
          ELEM_ACOUSTIC_IMPE_QUAD8R,
          ELEM_ACOUSTIC_IMPE_QUAD9,
          ELEM_ACOUSTIC_IMPE_TRIA6,
          ELEM_ACOUSTIC_IMPE_TRIA6H,
          ELEM_ACOUSTIC_IMPE_TRIA10,
          ELEM_ACOUSTIC_CONSTR,
          ELEM_ACOUSTIC_BND_ISO4,
          ELEM_ACOUSTIC_BND_QUAD8,
          ELEM_ACOUSTIC_BND_QUAD8R,
          ELEM_ACOUSTIC_BND_QUAD9,
          ELEM_ACOUSTIC_BND_TRIA6,
          ELEM_ACOUSTIC_BND_TRIA6H,
          ELEM_ACOUSTIC_BND_TRIA10,
          ELEM_FLUID_STRUCT_ISO4,
          ELEM_FLUID_STRUCT_QUAD8,
          ELEM_FLUID_STRUCT_QUAD8R,
          ELEM_FLUID_STRUCT_QUAD9,
          ELEM_FLUID_STRUCT_TRIA6,
          ELEM_FLUID_STRUCT_TRIA6H,
          ELEM_FLUID_STRUCT_TRIA10,
          ELEM_TYPE_COUNT,
          ELEM_TYPE_UNKNOWN = -1
     };

     struct TypeInfo {
          TypeId type;
          char name[16];
          octave_idx_type min_nodes, max_nodes;
          DofMap::ElementType dof_type;
     };

     static constexpr octave_idx_type iGetNumTypes() {
          return ELEM_TYPE_COUNT;
     }

     static const TypeInfo& GetType(octave_idx_type i) {
          FEM_ASSERT(i >= 0);
          FEM_ASSERT(i < iGetNumTypes());

          return rgElemTypes[i];
     }

private:
     static constexpr TypeInfo rgElemTypes[ELEM_TYPE_COUNT] = {
          {ElementTypes::ELEM_ISO8,                 "iso8",            8,  8, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_ISO20,                "iso20",          20, 20, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_ISO20R,               "iso20r",         20, 20, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_ISO27,                "iso27",          27, 27, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_PENTA15,              "penta15",        15, 15, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_TET10H,               "tet10h",         10, 10, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_TET10,                "tet10",          10, 10, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_TET20,                "tet20",          20, 20, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_BEAM2,                "beam2",           2,  2, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_RBE3,                 "rbe3",            2, -1, DofMap::ELEM_RBE3},
          {ElementTypes::ELEM_JOINT,                "joints",          1, -1, DofMap::ELEM_JOINT},
          {ElementTypes::ELEM_SPRING,               "springs",         1, -1, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_DASHPOT,              "dashpots",        1, -1, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_RIGID_BODY,           "bodies",          1,  1, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_SFNCON4,              "sfncon4",         1, -1, DofMap::ELEM_JOINT},
          {ElementTypes::ELEM_SFNCON6,              "sfncon6",         1, -1, DofMap::ELEM_JOINT},
          {ElementTypes::ELEM_SFNCON6H,             "sfncon6h",        1, -1, DofMap::ELEM_JOINT},
          {ElementTypes::ELEM_SFNCON8,              "sfncon8",         1, -1, DofMap::ELEM_JOINT},
          {ElementTypes::ELEM_SFNCON8R,             "sfncon8r",        1, -1, DofMap::ELEM_JOINT},
          {ElementTypes::ELEM_SFNCON9,              "sfncon9",         1, -1, DofMap::ELEM_JOINT},
          {ElementTypes::ELEM_SFNCON10,             "sfncon10",        1, -1, DofMap::ELEM_JOINT},
          {ElementTypes::ELEM_SFNCON4S,             "sfncon4s",        1, -1, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_SFNCON6S,             "sfncon6s",        1, -1, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_SFNCON6HS,            "sfncon6hs",       1, -1, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_SFNCON8S,             "sfncon8s",        1, -1, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_SFNCON8RS,            "sfncon8rs",       1, -1, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_SFNCON9S,             "sfncon9s",        1, -1, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_SFNCON10S,            "sfncon10s",       1, -1, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_PRESSURE_ISO4,        "iso4",            4,  4, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_PRESSURE_QUAD8,       "quad8",           8,  8, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_PRESSURE_QUAD8R,      "quad8r",          8,  8, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_PRESSURE_QUAD9,       "quad9",           9,  9, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_PRESSURE_TRIA6,       "tria6",           6,  6, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_PRESSURE_TRIA6H,      "tria6h",          6,  6, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_PRESSURE_TRIA10,      "tria10",         10, 10, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_STRUCT_FORCE,         "force",           1, -1, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_THERM_CONV_ISO4,      "iso4",            4,  4, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_THERM_CONV_QUAD8,     "quad8",           8,  8, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_THERM_CONV_QUAD8R,    "quad8r",          8,  8, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_THERM_CONV_QUAD9,     "quad9",           9,  9, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_THERM_CONV_TRIA6,     "tria6",           6,  6, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_THERM_CONV_TRIA6H,    "tria6h",          6,  6, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_THERM_CONV_TRIA10,    "tria10",         10, 10, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_THERM_CONSTR,         "thermal_constr",  1, -1, DofMap::ELEM_JOINT},
          {ElementTypes::ELEM_HEAT_SOURCE_ISO4,     "iso4",            4,  4, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_HEAT_SOURCE_QUAD8,    "quad8",           8,  8, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_HEAT_SOURCE_QUAD8R,   "quad8r",          8,  8, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_HEAT_SOURCE_QUAD9,    "quad9",           9,  9, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_HEAT_SOURCE_TRIA6,    "tria6",           6,  6, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_HEAT_SOURCE_TRIA6H,   "tria6h",          6,  6, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_HEAT_SOURCE_TRIA10,   "tria10",         10, 10, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_PARTICLE_VEL_ISO4,    "iso4",            4,  4, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_PARTICLE_VEL_QUAD8,   "quad8",           8,  8, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_PARTICLE_VEL_QUAD8R,  "quad8r",          8,  8, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_PARTICLE_VEL_QUAD9,   "quad9",           9,  9, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_PARTICLE_VEL_TRIA6,   "tria6",           6,  6, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_PARTICLE_VEL_TRIA6H,  "tria6h",          6,  6, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_PARTICLE_VEL_TRIA10,  "tria10",         10, 10, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_ACOUSTIC_IMPE_ISO4,   "iso4",            4,  4, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_ACOUSTIC_IMPE_QUAD8,  "quad8",           8,  8, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_ACOUSTIC_IMPE_QUAD8R, "quad8r",          8,  8, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_ACOUSTIC_IMPE_QUAD9,  "quad9",           9,  9, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_ACOUSTIC_IMPE_TRIA6,  "tria6",           6,  6, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_ACOUSTIC_IMPE_TRIA6H, "tria6h",          6,  6, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_ACOUSTIC_IMPE_TRIA10, "tria10",         10, 10, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_ACOUSTIC_CONSTR,      "acoustic_constr", 1, -1, DofMap::ELEM_JOINT},
          {ElementTypes::ELEM_ACOUSTIC_BND_ISO4,    "iso4",            4,  4, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_ACOUSTIC_BND_QUAD8,   "quad8",           8,  8, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_ACOUSTIC_BND_QUAD8R,  "quad8r",          8,  8, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_ACOUSTIC_BND_QUAD9,   "quad9",           9,  9, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_ACOUSTIC_BND_TRIA6,   "tria6",           6,  6, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_ACOUSTIC_BND_TRIA6H,  "tria6h",          6,  6, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_ACOUSTIC_BND_TRIA10,  "tria10",         10, 10, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_FLUID_STRUCT_ISO4,    "iso4",            4,  4, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_FLUID_STRUCT_QUAD8,   "quad8",           8,  8, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_FLUID_STRUCT_QUAD8R,  "quad8r",          8,  8, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_FLUID_STRUCT_QUAD9,   "quad9",           9,  9, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_FLUID_STRUCT_TRIA6,   "tria6",           6,  6, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_FLUID_STRUCT_TRIA6H,  "tria6h",          6,  6, DofMap::ELEM_NODOF},
          {ElementTypes::ELEM_FLUID_STRUCT_TRIA10,  "tria10",         10, 10, DofMap::ELEM_NODOF}
     };
};

constexpr ElementTypes::TypeInfo ElementTypes::rgElemTypes[ELEM_TYPE_COUNT];

struct PerfectlyMatchedLayer {
     explicit PerfectlyMatchedLayer(const octave_scalar_map& elements)
          : oData{GetPMLData(elements)} {
     }

     static octave_scalar_map GetPMLData(const octave_scalar_map& elements) {
          const auto iter_PML = elements.seek("perfectly_matched_layers");

          if (iter_PML == elements.end()) {
               return octave_scalar_map{};
          }

          const octave_value ov_PML = elements.contents(iter_PML);

          if (!(ov_PML.isstruct() && ov_PML.numel() == 1)) {
               throw std::runtime_error("PML domain: mesh.elements.perfectly_matched_layers must be a scalar struct");
          }

          return ov_PML.scalar_map_value();
     }

     const octave_scalar_map oData;

     struct ElemF {
          ComplexNDArray f;
          NDArray e1, e2;
     };

     array<ElemF, ElementTypes::ELEM_TYPE_COUNT> rgElem;
};

class Element
{
     enum MatTypeBits: unsigned {
          MAT_ID_SHIFT = 8u,
          MAT_ID_MASK = 0xFF00u
     };

public:
     enum MatSymmetryFlag: unsigned {
          MAT_SYM_FULL = 0x00u,
          MAT_SYM_UPPER = 0x01u,
          MAT_SYM_LOWER = 0x02u,
          MAT_SYM_DIAG = 0x03u,
          MAT_SYM_MASK = 0x03u,
          MAT_TYPE_MATRIX = 0x04u,
          MAT_TYPE_VECTOR = 0x08u,
          MAT_TYPE_SCALAR = 0x12u,
          MAT_TYPE_ARRAY = 0x16u,
          MAT_TYPE_MASK = 0x1cu,
          MAT_UPDATE_INFO_ALWAYS = 0x20u,
          MAT_COLL_PNT_OUTPUT = 0x40u
     };

     static_assert(((MAT_SYM_MASK | MAT_TYPE_MASK | MAT_UPDATE_INFO_ALWAYS) & MAT_ID_MASK) == 0u);
     static_assert(((MAT_SYM_MASK | MAT_TYPE_MASK | MAT_UPDATE_INFO_ALWAYS) & DofMap::DO_MASK) == 0u);
     static_assert((MAT_ID_MASK & DofMap::DO_MASK) == 0u);
     static_assert((MAT_SYM_MASK & MAT_TYPE_MASK) == 0u);
     static_assert((MAT_SYM_MASK & MAT_UPDATE_INFO_ALWAYS) == 0u);
     static_assert((MAT_TYPE_MASK & MAT_UPDATE_INFO_ALWAYS) == 0u);

     enum FemMatrixType: unsigned {
          MAT_UNKNOWN                    = 0u,
          MAT_STIFFNESS                  = ( 1u << MAT_ID_SHIFT) | MAT_TYPE_MATRIX | DofMap::DO_STRUCTURAL,
          MAT_STIFFNESS_SYM              =           MAT_STIFFNESS | MAT_SYM_UPPER,
          MAT_STIFFNESS_SYM_L            =           MAT_STIFFNESS | MAT_SYM_LOWER,
          MAT_MASS                       = ( 2u << MAT_ID_SHIFT) | MAT_TYPE_MATRIX | DofMap::DO_STRUCTURAL,
          MAT_MASS_SYM                   =           MAT_MASS | MAT_SYM_UPPER,
          MAT_MASS_SYM_L                 =           MAT_MASS | MAT_SYM_LOWER,
          MAT_MASS_LUMPED                =           MAT_MASS | MAT_SYM_DIAG,
          MAT_DAMPING                    = ( 3u << MAT_ID_SHIFT) | MAT_TYPE_MATRIX | DofMap::DO_STRUCTURAL,
          MAT_DAMPING_SYM                =           MAT_DAMPING | MAT_SYM_UPPER,
          MAT_DAMPING_SYM_L              =           MAT_DAMPING | MAT_SYM_LOWER,
          SCA_TOT_MASS                   = ( 4u << MAT_ID_SHIFT) | MAT_TYPE_SCALAR | DofMap::DO_STRUCTURAL,
          VEC_INERTIA_M1                 = ( 5u << MAT_ID_SHIFT) | MAT_TYPE_ARRAY | DofMap::DO_STRUCTURAL,
          MAT_INERTIA_J                  = ( 6u << MAT_ID_SHIFT) | MAT_TYPE_ARRAY | DofMap::DO_STRUCTURAL,
          MAT_INERTIA_INV3               = ( 7u << MAT_ID_SHIFT) | MAT_TYPE_ARRAY | DofMap::DO_STRUCTURAL,
          MAT_INERTIA_INV4               = ( 8u << MAT_ID_SHIFT) | MAT_TYPE_ARRAY | DofMap::DO_STRUCTURAL,
          MAT_INERTIA_INV5               = ( 9u << MAT_ID_SHIFT) | MAT_TYPE_ARRAY | DofMap::DO_STRUCTURAL,
          MAT_INERTIA_INV8               = (10u << MAT_ID_SHIFT) | MAT_TYPE_ARRAY | DofMap::DO_STRUCTURAL,
          MAT_INERTIA_INV9               = (11u << MAT_ID_SHIFT) | MAT_TYPE_ARRAY | DofMap::DO_STRUCTURAL,
          VEC_LOAD_CONSISTENT            = (12u << MAT_ID_SHIFT) | MAT_TYPE_VECTOR | DofMap::DO_STRUCTURAL,
          VEC_LOAD_LUMPED                = (13u << MAT_ID_SHIFT) | MAT_TYPE_VECTOR | DofMap::DO_STRUCTURAL,
          VEC_STRESS_CAUCH               = (14u << MAT_ID_SHIFT) | MAT_TYPE_ARRAY | DofMap::DO_STRUCTURAL,
          VEC_STRAIN_TOTAL               = (15u << MAT_ID_SHIFT) | MAT_TYPE_ARRAY | DofMap::DO_STRUCTURAL,
          SCA_STRESS_VMIS                = (16u << MAT_ID_SHIFT) | MAT_TYPE_ARRAY | DofMap::DO_STRUCTURAL,
          MAT_THERMAL_COND               = (17u << MAT_ID_SHIFT) | MAT_TYPE_MATRIX | DofMap::DO_THERMAL,
          MAT_HEAT_CAPACITY              = (18u << MAT_ID_SHIFT) | MAT_TYPE_MATRIX | DofMap::DO_THERMAL,
          VEC_LOAD_THERMAL               = (19u << MAT_ID_SHIFT) | MAT_TYPE_VECTOR | DofMap::DO_THERMAL,
          MAT_MASS_ACOUSTICS_RE          = (20u << MAT_ID_SHIFT) | MAT_TYPE_MATRIX | DofMap::DO_ACOUSTICS,
          MAT_MASS_ACOUSTICS_IM          = (21u << MAT_ID_SHIFT) | MAT_TYPE_MATRIX | DofMap::DO_ACOUSTICS,
          MAT_STIFFNESS_ACOUSTICS_RE     = (22u << MAT_ID_SHIFT) | MAT_TYPE_MATRIX | DofMap::DO_ACOUSTICS | MAT_UPDATE_INFO_ALWAYS,
          MAT_STIFFNESS_ACOUSTICS_IM     = (23u << MAT_ID_SHIFT) | MAT_TYPE_MATRIX | DofMap::DO_ACOUSTICS | MAT_UPDATE_INFO_ALWAYS,
          MAT_DAMPING_ACOUSTICS_RE       = (24u << MAT_ID_SHIFT) | MAT_TYPE_MATRIX | DofMap::DO_ACOUSTICS,
          MAT_DAMPING_ACOUSTICS_IM       = (25u << MAT_ID_SHIFT) | MAT_TYPE_MATRIX | DofMap::DO_ACOUSTICS,
          VEC_LOAD_ACOUSTICS             = (26u << MAT_ID_SHIFT) | MAT_TYPE_VECTOR | DofMap::DO_ACOUSTICS,
          VEC_PARTICLE_VELOCITY          = (27u << MAT_ID_SHIFT) | MAT_TYPE_ARRAY  | DofMap::DO_ACOUSTICS | DofMap::DO_FLUID_STRUCT,
          VEC_PARTICLE_VELOCITY_C        = (28u << MAT_ID_SHIFT) | MAT_TYPE_ARRAY  | DofMap::DO_ACOUSTICS | DofMap::DO_FLUID_STRUCT,
          SCA_ACOUSTIC_INTENSITY         = (29u << MAT_ID_SHIFT) | MAT_TYPE_ARRAY  | DofMap::DO_ACOUSTICS | DofMap::DO_FLUID_STRUCT,
          SCA_ACOUSTIC_INTENSITY_C       = (30u << MAT_ID_SHIFT) | MAT_TYPE_ARRAY  | DofMap::DO_ACOUSTICS | DofMap::DO_FLUID_STRUCT,
          VEC_SURFACE_NORMAL_VECTOR      = (31u << MAT_ID_SHIFT) | MAT_TYPE_ARRAY  | DofMap::DO_ACOUSTICS | DofMap::DO_FLUID_STRUCT,
          MAT_MASS_FLUID_STRUCT_RE       = (32u << MAT_ID_SHIFT) | MAT_TYPE_MATRIX | DofMap::DO_FLUID_STRUCT,
          MAT_MASS_FLUID_STRUCT_IM       = (33u << MAT_ID_SHIFT) | MAT_TYPE_MATRIX | DofMap::DO_FLUID_STRUCT,
          MAT_STIFFNESS_FLUID_STRUCT_RE  = (34u << MAT_ID_SHIFT) | MAT_TYPE_MATRIX | DofMap::DO_FLUID_STRUCT | MAT_UPDATE_INFO_ALWAYS,
          MAT_STIFFNESS_FLUID_STRUCT_IM  = (35u << MAT_ID_SHIFT) | MAT_TYPE_MATRIX | DofMap::DO_FLUID_STRUCT | MAT_UPDATE_INFO_ALWAYS,
          MAT_DAMPING_FLUID_STRUCT_RE    = (36u << MAT_ID_SHIFT) | MAT_TYPE_MATRIX | DofMap::DO_FLUID_STRUCT,
          MAT_DAMPING_FLUID_STRUCT_IM    = (37u << MAT_ID_SHIFT) | MAT_TYPE_MATRIX | DofMap::DO_FLUID_STRUCT,
          VEC_LOAD_FLUID_STRUCT          = (38u << MAT_ID_SHIFT) | MAT_TYPE_MATRIX | DofMap::DO_FLUID_STRUCT,
          MAT_STIFFNESS_TAU0             = (39u << MAT_ID_SHIFT) | MAT_TYPE_MATRIX | DofMap::DO_STRUCTURAL,
          MAT_STIFFNESS_OMEGA            = (40u << MAT_ID_SHIFT) | MAT_TYPE_MATRIX | DofMap::DO_STRUCTURAL,
          MAT_STIFFNESS_OMEGA_DOT        = (41u << MAT_ID_SHIFT) | MAT_TYPE_MATRIX | DofMap::DO_STRUCTURAL,
          MAT_DAMPING_OMEGA              = (42u << MAT_ID_SHIFT) | MAT_TYPE_MATRIX | DofMap::DO_STRUCTURAL,
          VEC_SURFACE_AREA               = (43u << MAT_ID_SHIFT) | MAT_TYPE_VECTOR | DofMap::DO_STRUCTURAL,
          MAT_STIFFNESS_IM               = (44u << MAT_ID_SHIFT) | MAT_TYPE_MATRIX | DofMap::DO_STRUCTURAL,
          VEC_COLL_MASS                  = MAT_MASS | MAT_COLL_PNT_OUTPUT,
          VEC_COLL_STIFFNESS             = MAT_STIFFNESS | MAT_COLL_PNT_OUTPUT,
          VEC_COLL_HEAT_CAPACITY         = MAT_HEAT_CAPACITY | MAT_COLL_PNT_OUTPUT,
          VEC_COLL_THERMAL_COND          = MAT_THERMAL_COND | MAT_COLL_PNT_OUTPUT,
          VEC_COLL_MASS_ACOUSTICS        = MAT_MASS_ACOUSTICS_RE | MAT_COLL_PNT_OUTPUT,
          VEC_COLL_STIFF_ACOUSTICS       = MAT_STIFFNESS_ACOUSTICS_RE | MAT_COLL_PNT_OUTPUT,
          VEC_COLL_MASS_FLUID_STRUCT     = MAT_MASS_FLUID_STRUCT_RE | MAT_COLL_PNT_OUTPUT,
          VEC_COLL_STIFF_FLUID_STRUCT    = MAT_STIFFNESS_FLUID_STRUCT_RE | MAT_COLL_PNT_OUTPUT
     };

     static constexpr unsigned MAT_TYPE_COUNT = 43u;

     static unsigned GetMatTypeIndex(FemMatrixType eMatType) {
          unsigned i = ((eMatType & MAT_ID_MASK) >> MAT_ID_SHIFT) - 1u;

          FEM_ASSERT(i < MAT_TYPE_COUNT);

          return i;
     }

     static constexpr FemMatrixType ComplexMatrixType(FemMatrixType eMatType) {
          switch (eMatType) {
          case VEC_PARTICLE_VELOCITY:
               return VEC_PARTICLE_VELOCITY_C;
          case SCA_ACOUSTIC_INTENSITY:
               return SCA_ACOUSTIC_INTENSITY_C;
          default:
               throw std::logic_error("element: requested complex matrix type does not exist");
          }
     }

     Element(ElementTypes::TypeId eltype, octave_idx_type id, const Matrix& X, const Material* material, const int32NDArray& nodes)
          :eltype(eltype), id(id), X(X), material(material), nodes(nodes), bElemAssDone{false} {

          FEM_ASSERT(X.columns() == nodes.numel());
     }

     Element(const Element& oElem)
          :eltype(oElem.eltype), id(oElem.id), X(oElem.X), material(oElem.material), nodes(oElem.nodes), bElemAssDone{false} {
     }

     virtual ~Element() {
     }

     virtual void Assemble(MatrixAss& mat, MeshInfo& info, const DofMap& dof, FemMatrixType eMatType) const=0;

     virtual void PostProcElem(FemMatrixType eMatType, PostProcData& oSolution) const {
     }

     virtual octave_idx_type iGetWorkSpaceSize(FemMatrixType eMatType) const=0;

     virtual double dGetMass() const {
          return 0.;
     }

     static constexpr bool bNeedMatrixInfo(FemMatrixType eMatType) {
          return false;
     }

     virtual void Extract(octave_idx_type& idx, octave_map& sElem) const {
     }

     virtual octave_idx_type iGetNumCollocPoints(FemMatrixType eMatType) const {
          return 0;
     }

     bool bSetElemAssDone() const {
          bool bExpected = false;
          return bElemAssDone.compare_exchange_strong(bExpected, true);
     }

     void ResetElemAssDone() const {
          bElemAssDone = false;
     }
protected:
     const ElementTypes::TypeId eltype;
     const octave_idx_type id;
     const Matrix X;
     const Material* material;
     const int32NDArray nodes;
     mutable std::atomic<bool> bElemAssDone;
};

class PostProcData {
public:
     enum FieldTypeReal: unsigned {
          VEC_NO_STRUCT_DISPLACEMENT_RE,     // nodal displacement vector field
          SCA_NO_ACOUSTIC_PART_VEL_POT_RE,   // nodal particle velocity potential scalar field
          SCA_NO_ACOUSTIC_PART_VEL_POT_P_RE, // first derivative of nodal particle velocity potential scalar field
          VEC_NO_ACOUSTIC_PART_VEL_RE,       // nodal particle velocity vector field
          SCA_NO_ACOUSTIC_PRESSURE_RE,       // nodal acoustic pressure scalar field
          VEC_EL_ACOUSTIC_PART_VEL_RE,       // particle velocity scalar field at elements
          SCA_EL_ACOUSTIC_PART_VEL_NORM_RE,  // surface normal particle velocity at elements
          SCA_EL_ACOUSTIC_INTENSITY_RE,      // acoustic intensity at elements
          SCA_EL_ACOUSTIC_SOUND_POWER_RE,    // acoustic sound power at elements
          VEC_EL_STRUCT_STRESS_CAUCH_RE,     // stress tensor at elements
          VEC_EL_STRUCT_STRAIN_TOTAL_RE,     // strain tensor at elements
          SCA_NO_THERMAL_TEMPERATURE_RE,     // scalar temperature field
          VEC_G_STRUCT_INERTIA_M1_RE,        // center of mass times total mass
          MAT_G_STRUCT_INERTIA_J_RE,         // total moment of inertia
          MAT_G_STRUCT_INERTIA_INV3_RE,      // MBDyn's modal invariant 3
          MAT_G_STRUCT_INERTIA_INV4_RE,      // MBDyn's modal invariant 4
          MAT_G_STRUCT_INERTIA_INV5_RE,      // MBDyn's modal invariant 5
          MAT_G_STRUCT_INERTIA_INV8_RE,      // MBDyn's modal invariant 8
          MAT_G_STRUCT_INERTIA_INV9_RE,      // MBDyn's modal invariant 9
          VEC_EL_SURFACE_NORMAL_VECTOR_RE,   // surface normal vector at elements
          VEC_EL_COLLOC_POINTS_RE,           // Cartesion coordinates of collocation points
          VEC_EL_SURFACE_AREA_RE,            // Surface area needed as weight factors for RBE3 elements
          FIELD_COUNT_RE,
          UNKNOWN_FIELD_RE = ~0u
     };

     enum FieldTypeComplex: unsigned {
          VEC_NO_STRUCT_DISPLACEMENT_C,
          SCA_NO_ACOUSTIC_PART_VEL_POT_C,
          SCA_NO_ACOUSTIC_PART_VEL_POT_P_C,
          VEC_NO_ACOUSTIC_PART_VEL_C,
          SCA_NO_ACOUSTIC_PRESSURE_C,
          VEC_EL_ACOUSTIC_PART_VEL_C,
          SCA_EL_ACOUSTIC_PART_VEL_NORM_C,
          SCA_EL_ACOUSTIC_INTENSITY_C,
          FIELD_COUNT_C,
          UNKNOWN_FIELD_C = ~0u
     };

     PostProcData(const DofMap& oDofMap, const octave_scalar_map& sol)
          :eDomain(oDofMap.GetDomain()), iNumSteps(0) {

          enum DataType {
               DT_REAL,
               DT_COMPLEX
          };

          struct SolField {
               unsigned id;
               DataType type;
               octave_idx_type cols;
               char name[6];
          };

          static constexpr SolField structFields[] = {
               {VEC_NO_STRUCT_DISPLACEMENT_RE, DT_REAL, 6, "def"},
               {VEC_NO_STRUCT_DISPLACEMENT_C, DT_COMPLEX, 6, "def"}
          };

          static constexpr SolField thermalFields[] = {
               {SCA_NO_THERMAL_TEMPERATURE_RE, DT_REAL, 1, "theta"}
          };

          static constexpr SolField acousticFields[] = {
               {SCA_NO_ACOUSTIC_PART_VEL_POT_RE, DT_REAL, 1, "Phi"},
               {SCA_NO_ACOUSTIC_PART_VEL_POT_C, DT_COMPLEX, 1, "Phi"},
               {SCA_NO_ACOUSTIC_PART_VEL_POT_P_RE, DT_REAL, 1, "PhiP"},
               {SCA_NO_ACOUSTIC_PART_VEL_POT_P_C, DT_COMPLEX, 1, "PhiP"}
          };

          static constexpr SolField fluidStructFields[] = {
               {VEC_NO_STRUCT_DISPLACEMENT_RE, DT_REAL, 6, "def"},
               {VEC_NO_STRUCT_DISPLACEMENT_C, DT_COMPLEX, 6, "def"},
               {SCA_NO_ACOUSTIC_PART_VEL_POT_RE, DT_REAL, 1, "Phi"},
               {SCA_NO_ACOUSTIC_PART_VEL_POT_C, DT_COMPLEX, 1, "Phi"},
               {SCA_NO_ACOUSTIC_PART_VEL_POT_P_RE, DT_REAL, 1, "PhiP"},
               {SCA_NO_ACOUSTIC_PART_VEL_POT_P_C, DT_COMPLEX, 1, "PhiP"}
          };

          struct DomainField {
               DofMap::DomainType domain;
               const SolField* begin;
               const SolField* end;
          };

          static constexpr DomainField domainFields[] = {
               {DofMap::DO_STRUCTURAL,   std::begin(structFields),      std::end(structFields)},
               {DofMap::DO_THERMAL,      std::begin(thermalFields),     std::end(thermalFields)},
               {DofMap::DO_ACOUSTICS,    std::begin(acousticFields),    std::end(acousticFields)},
               {DofMap::DO_FLUID_STRUCT, std::begin(fluidStructFields), std::end(fluidStructFields)}
          };

          constexpr octave_idx_type numDomainFields = sizeof(domainFields) / sizeof(domainFields[0]);

          octave_idx_type iDomain = -1;

          for (octave_idx_type i = 0; i < numDomainFields; ++i) {
               if (domainFields[i].domain == eDomain) {
                    iDomain = i;
                    break;
               }
          }

          if (iDomain < 0) {
               throw std::logic_error("post proc data: unknown domain for post processing");
          }

          const DomainField& oDomain = domainFields[iDomain];

          for (const SolField* pSol = oDomain.begin; pSol != oDomain.end; ++pSol) {
               const auto iterSol = sol.seek(pSol->name);

               if (iterSol == sol.end()) {
                    continue;
               }

               const octave_value ovSol = sol.contents(iterSol);

               if (!(ovSol.is_matrix_type() && (ovSol.isreal() || ovSol.iscomplex()))) {
                    throw std::runtime_error("post proc data: field sol."s + pSol->name + " must be an real or complex array in argument sol");
               }

               dim_vector dimSol;

               if (ovSol.isreal() && pSol->type == DT_REAL) {
                    const NDArray solArrayRe = ovSol.array_value();
                    dimSol = solArrayRe.dims();
                    SetField(static_cast<FieldTypeReal>(pSol->id), ElementTypes::ELEM_TYPE_UNKNOWN, solArrayRe);
               } else if (ovSol.iscomplex() && pSol->type == DT_COMPLEX) {
                    const ComplexNDArray solArrayC = ovSol.complex_array_value();
                    dimSol = solArrayC.dims();
                    SetField(static_cast<FieldTypeComplex>(pSol->id), ElementTypes::ELEM_TYPE_UNKNOWN, solArrayC);
               } else {
                    continue;
               }

               octave_idx_type iNumCols, iNumStepsCurr;

               if (pSol->cols == 1) {
                    iNumCols = pSol->cols;
                    iNumStepsCurr = dimSol(1);
               } else {
                    iNumCols = dimSol(1);
                    iNumStepsCurr = dimSol.ndims() > 2 ? dimSol(2) : 1;
               }

               if (iNumCols != pSol->cols) {
                    throw std::runtime_error("post proc data: columns of sol."s + pSol->name + " is not valid");
               }

               if (iNumSteps > 0 && iNumStepsCurr != iNumSteps) {
                    throw std::runtime_error("post proc data: number of load steps in sol."s + pSol->name + " is not consistent within sol");
               }

               iNumSteps = iNumStepsCurr;
          }
     }

     NDArray& GetField(FieldTypeReal eFieldType, ElementTypes::TypeId eElemType) {
          auto iter = nodalFieldsReal.find(Key<FieldTypeReal>{eFieldType, eElemType});

          if (iter == nodalFieldsReal.end()) {
               throw std::runtime_error("post proc data: real field "s + szFieldName[eFieldType] + " not found in argument sol");
          }

          return iter->second;
     }

     ComplexNDArray& GetField(FieldTypeComplex eFieldType, ElementTypes::TypeId eElemType) {
          auto iter = nodalFieldsComplex.find(Key<FieldTypeComplex>{eFieldType, eElemType});

          if (iter == nodalFieldsComplex.end()) {
               throw std::runtime_error("post proc data: complex field "s + szFieldName[eFieldType] + " not found in argument sol");
          }

          return iter->second;
     }

     void SetField(FieldTypeReal eFieldType, ElementTypes::TypeId eElemType, const NDArray& fieldData) {
          nodalFieldsReal[Key<FieldTypeReal>{eFieldType, eElemType}] = fieldData;
     }

     void SetField(FieldTypeComplex eFieldType, ElementTypes::TypeId eElemType, const ComplexNDArray& fieldData) {
          nodalFieldsComplex[Key<FieldTypeComplex>{eFieldType, eElemType}] = fieldData;
     }

     DofMap::DomainType GetDomain() const { return eDomain; }

     octave_idx_type GetNumSteps() const { return iNumSteps; }

     static bool bIsNodalField(FieldTypeReal fieldType) {
          switch (fieldType) {
          case VEC_NO_STRUCT_DISPLACEMENT_RE:
          case SCA_NO_ACOUSTIC_PART_VEL_POT_RE:
          case SCA_NO_ACOUSTIC_PART_VEL_POT_P_RE:
          case VEC_NO_ACOUSTIC_PART_VEL_RE:
          case SCA_NO_ACOUSTIC_PRESSURE_RE:
          case SCA_NO_THERMAL_TEMPERATURE_RE:
          case VEC_G_STRUCT_INERTIA_M1_RE:
          case MAT_G_STRUCT_INERTIA_J_RE:
          case MAT_G_STRUCT_INERTIA_INV3_RE:
          case MAT_G_STRUCT_INERTIA_INV4_RE:
          case MAT_G_STRUCT_INERTIA_INV5_RE:
          case MAT_G_STRUCT_INERTIA_INV8_RE:
          case MAT_G_STRUCT_INERTIA_INV9_RE:
               return true;
          default:
               return false;
          }
     }

     static bool bIsNodalField(FieldTypeComplex fieldType) {
          switch (fieldType) {
          case VEC_NO_STRUCT_DISPLACEMENT_C:
          case SCA_NO_ACOUSTIC_PART_VEL_POT_C:
          case SCA_NO_ACOUSTIC_PART_VEL_POT_P_C:
          case VEC_NO_ACOUSTIC_PART_VEL_C:
          case SCA_NO_ACOUSTIC_PRESSURE_C:
               return true;
          default:
               return false;
          }
     }

     static constexpr FieldTypeComplex ComplexFieldType(FieldTypeReal fieldType) {
          const unsigned index = fieldType;

          if (!(index < FIELD_COUNT_C)) {
               throw std::logic_error("post proc data: requested complex postprocessing field does not exist");
          }

          return static_cast<FieldTypeComplex>(fieldType);
     }

private:
     template <typename FieldType>
     struct Key {
          Key(FieldType fieldType, ElementTypes::TypeId elemType)
               :field(fieldType), element(elemType) {

               if (PostProcData::bIsNodalField(field)) {
                    element = ElementTypes::ELEM_TYPE_UNKNOWN;
               }
          }

          bool operator==(const Key& key) const {
               return field == key.field && element == key.element;
          }

          FieldType field;
          ElementTypes::TypeId element;
     };

     template <typename FieldType>
     struct Hash {
          std::size_t operator()(const Key<FieldType>& key) const noexcept {
               std::size_t h1 = std::hash<unsigned>{}(key.field);
               std::size_t h2 = std::hash<unsigned>{}(key.element);
               return h1 ^ (h2 << 1);
          }
     };

     const DofMap::DomainType eDomain;
     std::unordered_map<Key<FieldTypeReal>, NDArray, Hash<FieldTypeReal> > nodalFieldsReal;
     std::unordered_map<Key<FieldTypeComplex>, ComplexNDArray, Hash<FieldTypeComplex> > nodalFieldsComplex;
     octave_idx_type iNumSteps;
     static constexpr char szFieldName[FIELD_COUNT_RE][8] = {"def",
                                                             "Phi",
                                                             "PhiP",
                                                             "v",
                                                             "p",
                                                             "v",
                                                             "vn",
                                                             "tau",
                                                             "epsilon",
                                                             "theta",
                                                             "M1",
                                                             "J",
                                                             "Inv3",
                                                             "Inv4",
                                                             "Inv5",
                                                             "Inv8",
                                                             "Inv9"};
};

template <typename T>
struct PostProcFieldHelper;

template <>
struct PostProcFieldHelper<double> {
     static constexpr PostProcData::FieldTypeReal
     ConvertFieldType(PostProcData::FieldTypeReal eFieldTypeReal) {
          return eFieldTypeReal;
     }

     static constexpr Element::FemMatrixType
     ConvertMatrixType(Element::FemMatrixType eMatType) {
          return eMatType;
     }
};

template <>
struct PostProcFieldHelper<std::complex<double> > {
     static constexpr PostProcData::FieldTypeComplex
     ConvertFieldType(PostProcData::FieldTypeReal eFieldTypeReal) {
          return PostProcData::ComplexFieldType(eFieldTypeReal);
     }

     static constexpr Element::FemMatrixType
     ConvertMatrixType(Element::FemMatrixType eMatType) {
          return Element::ComplexMatrixType(eMatType);
     }
};

constexpr char PostProcData::szFieldName[FIELD_COUNT_RE][8];

class MatrixAss {
public:
     struct MatrixInfo {
          double beta = 1.; // scale factor for constraint equations
          bool updated = false;
     };

     explicit MatrixAss(octave_idx_type max_nnz)
          :eMatType(Element::MAT_UNKNOWN),
           nnz(0),
           ridx(dim_vector(max_nnz, 1), -1),
           cidx(dim_vector(max_nnz, 1), -1),
           data(max_nnz, 0.) {
     }

     bool bHaveMatrixInfo(Element::FemMatrixType eMatType) const {
          unsigned idx = Element::GetMatTypeIndex(eMatType);

          return info[idx].updated;
     }

     const MatrixInfo& GetMatrixInfo(Element::FemMatrixType eMatType) const {
          unsigned idx = Element::GetMatTypeIndex(eMatType);

          if (!info[idx].updated) {
               throw std::runtime_error("fem_ass_matrix: invalid order of values in argument matrix_type - e.g. the right hand side vector must be assembled always after the stiffness matrix");
          }

          return info[idx];
     }

     void UpdateMatrixInfo(const DofMap& oDofMap) {
          const unsigned iMatType = Element::GetMatTypeIndex(eMatType);

          if (info[iMatType].updated) {
               return;
          }

          ColumnVector diagA(oDofMap.iGetNumDof(), 0.);

          for (octave_idx_type i = 0; i < nnz; ++i) {
               if (ridx.xelem(i) == cidx.xelem(i)) {
                    diagA.xelem(ridx.xelem(i)) += data.xelem(i);
               }
          }

          constexpr double INIT_MIN = std::numeric_limits<double>::max();
          constexpr double INIT_MAX = -std::numeric_limits<double>::max();

          double minA = INIT_MIN;
          double maxA = INIT_MAX;

          for (octave_idx_type i = 0; i < diagA.numel(); ++i) {
               const double absA = fabs(diagA.xelem(i));

               if (absA > maxA) {
                    maxA = absA;
               }

               if (absA < minA) {
                    minA = absA;
               }
          }

          // According to Code_Aster r3.03.08
          info[iMatType].beta = (minA > 0. || maxA > 0.) ? 0.5 * (minA + maxA) : 1.;
          info[iMatType].updated = true;
     }

     void Insert(double d, octave_idx_type r, octave_idx_type c) {
          if (bNeedToInsertElem(r, c)) {
               InsertRaw(d, r, c);
          }
     }

     void Insert(const Matrix& Ke, const Array<octave_idx_type>& r, const Array<octave_idx_type>& c) {
          const octave_idx_type nrows = Ke.rows();
          const octave_idx_type ncols = Ke.columns();

          for (octave_idx_type j = 0; j < ncols; ++j) {
               for (octave_idx_type i = 0; i < nrows; ++i) {
                    Insert(Ke.xelem(i + nrows * j), r.xelem(i), c.xelem(j));
               }
          }
     }

     void Finish() {
          // Do not resize the workspace here because it could be reused for other matrices!
     }

     SparseMatrix Assemble(const DofMap& oDofMap, octave_idx_type iNumLoads) const {
          FEM_ASSERT(eMatType != Element::MAT_UNKNOWN);

          const octave_idx_type iNumRows = oDofMap.iGetNumDof();
          octave_idx_type iNumCols = -1;

          switch (eMatType) {
          case Element::VEC_LOAD_CONSISTENT:
          case Element::VEC_LOAD_LUMPED:
          case Element::VEC_LOAD_THERMAL:
          case Element::VEC_LOAD_ACOUSTICS:
          case Element::VEC_LOAD_FLUID_STRUCT:
               iNumCols = iNumLoads;
               break;
#if DEBUG > 0
          case Element::MAT_UNKNOWN:
               FEM_ASSERT(0);
               break;
#endif
          default:
               iNumCols = iNumRows;
          }

          return SparseMatrix(data.linear_slice(0, nnz),
                              ridx.linear_slice(0, nnz),
                              cidx.linear_slice(0, nnz),
                              iNumRows,
                              iNumCols);
     }

     void Reset(Element::FemMatrixType eMatTypeCurr) {
          eMatType = eMatTypeCurr;
          nnz = 0;
     }

private:
     bool bNeedToInsertElem(octave_idx_type r, octave_idx_type c) const {
          if (!(r > 0 && c > 0)) {
               return false;
          }

          switch (eMatType & Element::MAT_SYM_MASK) {
          case Element::MAT_SYM_UPPER:
               if (c < r) {
                    return false;
               }
               break;
          case Element::MAT_SYM_LOWER:
               if (c > r) {
                    return false;
               }
               break;
          case Element::MAT_SYM_DIAG:
               if (c != r) {
                    return false;
               }
               break;
          default:
               break;
          }

          return true;
     }

     void InsertRaw(double d, octave_idx_type r, octave_idx_type c) {
          FEM_ASSERT(bNeedToInsertElem(r, c));

          const octave_idx_type current = nnz++;

          if (data.rows() <= current) {
               throw std::runtime_error("fem_ass_matrix: allocated workspace size exceeded");
          }

          data.xelem(current) = d;
          ridx.xelem(current) = r - 1;
          cidx.xelem(current) = c - 1;
     }

     Element::FemMatrixType eMatType;
     std::atomic<octave_idx_type> nnz;
     Array<octave_idx_type> ridx, cidx;
     ColumnVector data;
     std::array<MatrixInfo, Element::MAT_TYPE_COUNT> info;
};

class ElemJoint: public Element
{
public:
     ElemJoint(ElementTypes::TypeId eltype, octave_idx_type id, const Matrix& X, const Material* material, const int32NDArray& nodes, const Matrix& C, const Matrix& U, DofMap::DomainType eDomain, double dScale)
          :Element(eltype, id, X, material, nodes),
           C(C),
           U(U),
           iNumNodeDof(eltype == ElementTypes::ELEM_JOINT ? 6 : 1),
           eNodalDofType(GetNodalDofType(eltype)),
           dScale(dScale) {
          FEM_ASSERT(C.columns() == nodes.numel() * iNumNodeDof);
          FEM_ASSERT(C.rows() <= C.columns());
          FEM_ASSERT(C.rows() >= 1);
          FEM_ASSERT(U.rows() == C.rows());
          FEM_ASSERT(X.rows() == 6);
     }

     virtual void Assemble(MatrixAss& mat, MeshInfo& info, const DofMap& dof, const FemMatrixType eMatType) const override {
          double beta;
          const octave_idx_type Crows = C.rows();
          const octave_idx_type Ccols = C.columns();

          switch (eMatType) {
          case MAT_STIFFNESS:
          case MAT_STIFFNESS_SYM:
          case MAT_STIFFNESS_SYM_L:
          case VEC_LOAD_CONSISTENT:
          case VEC_LOAD_LUMPED:
               beta = mat.GetMatrixInfo(MAT_STIFFNESS).beta;
               break;

          case MAT_THERMAL_COND:
          case VEC_LOAD_THERMAL:
               beta = mat.GetMatrixInfo(MAT_THERMAL_COND).beta;
               break;

          case MAT_DAMPING_ACOUSTICS_RE:
          case VEC_LOAD_ACOUSTICS:
               beta = mat.GetMatrixInfo(MAT_STIFFNESS_ACOUSTICS_RE).beta;
               break;

          case MAT_STIFFNESS_FLUID_STRUCT_RE:
          case MAT_DAMPING_FLUID_STRUCT_RE:
          case VEC_LOAD_FLUID_STRUCT:
               beta = mat.GetMatrixInfo(MAT_STIFFNESS_FLUID_STRUCT_RE).beta;
               break;

          default:
               beta = 1.;
          }

          beta *= dScale;

          switch (eMatType) {
          case MAT_STIFFNESS:
          case MAT_STIFFNESS_SYM:
          case MAT_STIFFNESS_SYM_L:
          case MAT_THERMAL_COND:
          case MAT_DAMPING_ACOUSTICS_RE:
          case MAT_DAMPING_FLUID_STRUCT_RE:
          case MAT_STIFFNESS_FLUID_STRUCT_RE:
          {
               if (eMatType == MAT_STIFFNESS_FLUID_STRUCT_RE && eNodalDofType != DofMap::NDOF_DISPLACEMENT) {
                    return;
               }

               if (eMatType == MAT_DAMPING_FLUID_STRUCT_RE && eNodalDofType != DofMap::NDOF_VELOCITY_POT) {
                    return;
               }

               Array<octave_idx_type> ndofidx(dim_vector(nodes.numel() * iNumNodeDof, 1), -1);

               for (octave_idx_type inode = 0; inode < nodes.numel(); ++inode) {
                    for (octave_idx_type idof = 0; idof < iNumNodeDof; ++idof) {
                         ndofidx.xelem(inode * iNumNodeDof + idof) = dof.GetNodeDofIndex(nodes.xelem(inode).value() - 1, eNodalDofType, idof);
                    }
               }

               Array<octave_idx_type> edofidx(dim_vector(Crows, 1));

               for (octave_idx_type idof = 0; idof < edofidx.numel(); ++idof) {
                    edofidx.xelem(idof) = dof.GetElemDofIndex(DofMap::ELEM_JOINT, id - 1, idof);
               }

               for (octave_idx_type j = 0; j < Ccols; ++j) {
                    for (octave_idx_type i = 0; i < Crows; ++i) {
                         const double Cij = beta * C.xelem(i + Crows * j);
                         mat.Insert(Cij, ndofidx.xelem(j), edofidx.xelem(i));
                         mat.Insert(Cij, edofidx.xelem(i), ndofidx.xelem(j));
                    }
               }
          } break;
          case VEC_LOAD_CONSISTENT:
          case VEC_LOAD_LUMPED:
          case VEC_LOAD_THERMAL:
          case VEC_LOAD_ACOUSTICS:
          case VEC_LOAD_FLUID_STRUCT: {
               Array<octave_idx_type> edofidx(dim_vector(Crows, 1));

               for (octave_idx_type idof = 0; idof < Crows; ++idof) {
                    edofidx.xelem(idof) = dof.GetElemDofIndex(DofMap::ELEM_JOINT, id - 1, idof);
               }

               if (eNodalDofType == DofMap::NDOF_VELOCITY_POT) {
                    beta = -beta;
               }

               const octave_idx_type Ucols = U.columns();
               const octave_idx_type Urows = U.rows();

               for (octave_idx_type j = 0; j < Ucols; ++j) {
                    for (octave_idx_type i = 0; i < Urows; ++i) {
                         mat.Insert(beta * U.xelem(i + Urows * j), edofidx.xelem(i), j + 1);
                    }
               }
          } break;
          default:
               break;
          }
     }

     virtual octave_idx_type iGetWorkSpaceSize(FemMatrixType eMatType) const override {
          switch (eMatType) {
          case MAT_STIFFNESS:
          case MAT_STIFFNESS_SYM:
          case MAT_STIFFNESS_SYM_L:
          case MAT_STIFFNESS_FLUID_STRUCT_RE:
          case MAT_THERMAL_COND:
          case MAT_DAMPING_ACOUSTICS_RE:
          case MAT_DAMPING_FLUID_STRUCT_RE:
               return 2 * C.rows() * C.columns();
          case VEC_LOAD_CONSISTENT:
          case VEC_LOAD_LUMPED:
          case VEC_LOAD_THERMAL:
          case VEC_LOAD_ACOUSTICS:
               return U.rows() * U.columns();
          default:
               return 0;
          }
     }

     static constexpr bool bNeedMatrixInfo(Element::FemMatrixType eMatType) {
          switch (eMatType) {
          case MAT_STIFFNESS:
          case MAT_STIFFNESS_SYM:
          case MAT_STIFFNESS_SYM_L:
          case MAT_STIFFNESS_FLUID_STRUCT_RE:
          case MAT_THERMAL_COND:
          case MAT_DAMPING_ACOUSTICS_RE:
          case MAT_DAMPING_FLUID_STRUCT_RE:
          case VEC_LOAD_CONSISTENT:
          case VEC_LOAD_LUMPED:
          case VEC_LOAD_THERMAL:
          case VEC_LOAD_ACOUSTICS:
          case VEC_LOAD_FLUID_STRUCT:
               return true;
          default:
               return false;
          }
     }

     virtual void Extract(octave_idx_type& idx, octave_map& sElem) const override {
          FEM_ASSERT(sElem.numel() > idx);

          Cell& ovC = sElem.contents("C");
          Cell& ovNodes = sElem.contents("nodes");

          ovC(idx) = C;
          ovNodes(idx) = nodes.transpose();
          ++idx;
     }

     static octave_idx_type iGetNumDofNodeMax(ElementTypes::TypeId eltype) {
          switch (eltype) {
          case ElementTypes::ELEM_JOINT:
               return 6;
          default:
               return 1;
          }
     }

     static void AllocIntegrationRule(Element::FemMatrixType eMatType) {
     }
private:
     static DofMap::NodalDofType GetNodalDofType(const ElementTypes::TypeId eltype) {
          switch (eltype) {
          case ElementTypes::ELEM_JOINT:
               return DofMap::NDOF_DISPLACEMENT;

          case ElementTypes::ELEM_THERM_CONSTR:
               return DofMap::NDOF_TEMPERATURE;

          case ElementTypes::ELEM_ACOUSTIC_CONSTR:
               return DofMap::NDOF_VELOCITY_POT;

          default:
               throw std::logic_error("joint: nodal constraint type not supported");
          }
     }

     const Matrix C, U;
     const octave_idx_type iNumNodeDof;
     const DofMap::NodalDofType eNodalDofType;
     const double dScale;
};

template <typename ScalarType, typename MatType>
class ElemSpringDashpotBase: public Element
{
public:
     ElemSpringDashpotBase(ElementTypes::TypeId eltype, octave_idx_type id, const Matrix& X, const Material* material, const int32NDArray& nodes)
          :Element(eltype, id, X, material, nodes) {
          FEM_ASSERT(X.rows() == 6);
     }

     static constexpr bool bNeedMatrixInfo(Element::FemMatrixType eMatType) {
          return false;
     }

     static void AllocIntegrationRule(Element::FemMatrixType eMatType) {
     }

protected:
     void AssembleMatrixReal(MatrixAss& mat, const DofMap& dof, const MatType& A) const {
          AssembleMatrix<RealPart>(mat, dof, A);
     }

     void AssembleMatrixImag(MatrixAss& mat, const DofMap& dof, const MatType& A) const {
          AssembleMatrix<ImagPart>(mat, dof, A);
     }

     template <double ScalarFunctionRealImag(const ScalarType&)>
     void AssembleMatrix(MatrixAss& mat, const DofMap& dof, const MatType& A) const {
          constexpr octave_idx_type iNumNodeDof = 6;
          const octave_idx_type Arows = A.rows();
          const octave_idx_type Acols = A.columns();

          Array<octave_idx_type> ndofidx(dim_vector(nodes.numel() * iNumNodeDof, 1), -1);

          FEM_ASSERT(Arows == ndofidx.numel());
          FEM_ASSERT(Acols == ndofidx.numel());

          for (octave_idx_type inode = 0; inode < nodes.numel(); ++inode) {
               for (octave_idx_type idof = 0; idof < iNumNodeDof; ++idof) {
                    ndofidx.xelem(inode * iNumNodeDof + idof) = dof.GetNodeDofIndex(nodes.xelem(inode).value() - 1, DofMap::NDOF_DISPLACEMENT, idof);
               }
          }

          for (octave_idx_type j = 0; j < Acols; ++j) {
               for (octave_idx_type i = 0; i < Arows; ++i) {
                    const double Aij = ScalarFunctionRealImag(A.xelem(i + Arows * j));
                    mat.Insert(Aij, ndofidx.xelem(j), ndofidx.xelem(i));
               }
          }
     }

     static double RealPart(const ScalarType& z) {
          if constexpr(std::is_same<ScalarType, double>::value) {
               return z;
          } else {
               return std::real(z);
          }
     }

     static double ImagPart(const ScalarType& z) {
          return std::imag(z);
     }
};

template <typename ScalarType, typename MatType>
class ElemSpring: public ElemSpringDashpotBase<ScalarType, MatType>
{
     typedef ElemSpringDashpotBase<ScalarType, MatType> BaseType;
     using BaseType::AssembleMatrixReal;
     using BaseType::AssembleMatrixImag;
public:
     ElemSpring(ElementTypes::TypeId eltype, octave_idx_type id, const Matrix& X, const Material* material, const int32NDArray& nodes, const MatType& K)
          :BaseType(eltype, id, X, material, nodes), K(K) {
          FEM_ASSERT(K.rows() == nodes.numel() * 6);
          FEM_ASSERT(K.columns() == nodes.numel() * 6);
     }

     virtual void Extract(octave_idx_type& idx, octave_map& sElem) const override {
          FEM_ASSERT(sElem.numel() > idx);

          Cell& ovK = sElem.contents("K");
          Cell& ovNodes = sElem.contents("nodes");

          ovK(idx) = K;
          ovNodes(idx) = this->nodes.transpose();
          ++idx;
     }

     virtual void Assemble(MatrixAss& mat, MeshInfo& info, const DofMap& dof, const Element::FemMatrixType eMatType) const override {
          switch (eMatType) {
          case Element::MAT_STIFFNESS:
          case Element::MAT_STIFFNESS_SYM:
          case Element::MAT_STIFFNESS_SYM_L:
          case Element::MAT_STIFFNESS_FLUID_STRUCT_RE:
               AssembleMatrixReal(mat, dof, K);
               break;

          case Element::MAT_STIFFNESS_IM:
          case Element::MAT_STIFFNESS_FLUID_STRUCT_IM:
               if constexpr(std::is_same<ScalarType, std::complex<double>>::value) {
                   AssembleMatrixImag(mat, dof, K);
               }
               break;

          default:
               break;
          }
     }

     virtual octave_idx_type iGetWorkSpaceSize(Element::FemMatrixType eMatType) const override {
          switch (eMatType) {
          case Element::MAT_STIFFNESS:
          case Element::MAT_STIFFNESS_IM:
          case Element::MAT_STIFFNESS_SYM:
          case Element::MAT_STIFFNESS_SYM_L:
          case Element::MAT_STIFFNESS_FLUID_STRUCT_RE:
          case Element::MAT_STIFFNESS_FLUID_STRUCT_IM:
               return K.rows() * K.columns();
          default:
               return 0;
          }
     }
private:
     const MatType K;
};

template <typename ScalarType, typename MatType>
class ElemDashpot: public ElemSpringDashpotBase<ScalarType, MatType>
{
     typedef ElemSpringDashpotBase<ScalarType, MatType> BaseType;
     using BaseType::AssembleMatrixReal;
     using BaseType::AssembleMatrixImag;
public:
     ElemDashpot(ElementTypes::TypeId eltype, octave_idx_type id, const Matrix& X, const Material* material, const int32NDArray& nodes, const MatType& D)
          :BaseType(eltype, id, X, material, nodes), D(D) {
          FEM_ASSERT(D.rows() == nodes.numel() * 6);
          FEM_ASSERT(D.columns() == nodes.numel() * 6);
     }

     virtual void Extract(octave_idx_type& idx, octave_map& sElem) const override {
          FEM_ASSERT(sElem.numel() > idx);

          Cell& ovD = sElem.contents("D");
          Cell& ovNodes = sElem.contents("nodes");

          ovD(idx) = D;
          ovNodes(idx) = this->nodes.transpose();
          ++idx;
     }

     virtual void Assemble(MatrixAss& mat, MeshInfo& info, const DofMap& dof, const Element::FemMatrixType eMatType) const override {
          switch (eMatType) {
          case Element::MAT_DAMPING:
          case Element::MAT_DAMPING_SYM:
          case Element::MAT_DAMPING_SYM_L:
          case Element::MAT_DAMPING_FLUID_STRUCT_RE:
               AssembleMatrixReal(mat, dof, D);
               break;

               //case MAT_DAMPING_IM:
          case Element::MAT_DAMPING_FLUID_STRUCT_IM:
               if constexpr(std::is_same<ScalarType, std::complex<double>>::value) {
                   AssembleMatrixImag(mat, dof, D);
               }
               break;

          default:
               break;
          }
     }

     virtual octave_idx_type iGetWorkSpaceSize(Element::FemMatrixType eMatType) const override {
          switch (eMatType) {
          case Element::MAT_DAMPING:
               //case MAT_DAMPING_IM:
          case Element::MAT_DAMPING_SYM:
          case Element::MAT_DAMPING_SYM_L:
          case Element::MAT_DAMPING_FLUID_STRUCT_RE:
          case Element::MAT_DAMPING_FLUID_STRUCT_IM:
               return D.rows() * D.columns();
          default:
               return 0;
          }
     }

private:
     const MatType D;
};

template <typename ScalarType, typename MatType>
class ElemContact: public ElemSpring<ScalarType, MatType>
{
     typedef ElemSpring<ScalarType, MatType> BaseType;
public:
     ElemContact(ElementTypes::TypeId eltype, octave_idx_type id, const Matrix& X, const Material* material, const int32NDArray& nodes, const Matrix& R, const Matrix& Hf, const octave_value& fnk, const double h0, const double dA)
          :BaseType(eltype, id, X, material, nodes, StiffnessMatrix(R, Hf, fnk, h0, dA)) {
     }

     static MatType StiffnessMatrix(const Matrix& R, const Matrix& Hf, const octave_value& fnk, const double h0, const double dA) {
          octave_value_list args;

          args.append(R);
          args.append(h0);
          args.append(dA);

          octave_value_list ovk = octave::feval(fnk, args, 1);

          if (ovk.numel() != 1) {
               throw std::runtime_error("sfncon.k must return one output argument");
          }

          MatType km;

          if constexpr(std::is_same<ScalarType, double>::value) {
               km = ovk.xelem(0).matrix_value();
          } else {
               km = ovk.xelem(0).complex_matrix_value();
          }

          if (!(km.rows() == 3 && km.columns() == 3)) {
               throw std::runtime_error("sfncon.k must be a 3x3 matrix");
          }

          MatType kmR_T(3, 3);

          for (octave_idx_type j = 0; j < 3; ++j) {
               for (octave_idx_type i = 0; i < 3; ++i) {
                    ScalarType aij{};

                    for (octave_idx_type k = 0; k < 3; ++k) {
                         aij += km.xelem(i, k) * R.xelem(j, k);
                    }

                    kmR_T.xelem(i, j) = aij;
               }
          }

          MatType Khtau(3, 3);

          for (octave_idx_type j = 0; j < 3; ++j) {
               for (octave_idx_type i = 0; i < 3; ++i) {
                    ScalarType aij{};

                    for (octave_idx_type k = 0; k < 3; ++k) {
                         aij += R.xelem(i, k) * kmR_T.xelem(k, j);
                    }

                    Khtau.xelem(i, j) = aij * dA;
               }
          }

          const octave_idx_type M = Hf.columns();
          const octave_idx_type N = M / 3;

          FEM_ASSERT(M % 3 == 0);
          FEM_ASSERT(N > 0);
          FEM_ASSERT(N * 3 == M);
          FEM_ASSERT(Hf.rows() == 3);

          MatType Khtau_Hf(3, M);

          for (octave_idx_type j = 0; j < M; ++j) {
               for (octave_idx_type i = 0; i < 3; ++i) {
                    ScalarType aij{};

                    for (octave_idx_type k = 0; k < 3; ++k) {
                         aij += Khtau.xelem(i, k) * Hf.xelem(k, j);
                    }

                    Khtau_Hf.xelem(i, j) = aij;
               }
          }

          MatType KHtau((N + 1) * 6, (N + 1) * 6, 0.);

          for (octave_idx_type j = 0; j < 3; ++j) {
               for (octave_idx_type i = 0; i < 3; ++i) {
                    KHtau.xelem(i, j) = Khtau.xelem(i, j);
               }
          }

          for (octave_idx_type j = 0; j < N; ++j) {
               for (octave_idx_type jj = 0; jj < 3; ++jj) {
                    for (octave_idx_type i = 0; i < 3; ++i) {
                         KHtau.xelem(i, 6 * (j + 1) + jj) = KHtau.xelem(6 * (j + 1) + jj, i) = -Khtau_Hf.xelem(i, j * 3 + jj);
                    }
               }
          }

          for (octave_idx_type j = 0; j < N; ++j) {
               for (octave_idx_type jj = 0; jj < 3; ++jj) {
                    for (octave_idx_type i = 0; i < N; ++i) {
                         for (octave_idx_type ii = 0; ii < 3; ++ii) {
                              ScalarType aij{};

                              for (octave_idx_type k = 0; k < 3; ++k) {
                                   aij += Hf.xelem(k, i * 3 + ii) * Khtau_Hf.xelem(k, j * 3 + jj);
                              }

                              KHtau.xelem(6 * (i + 1) + ii, 6 * (j + 1) + jj) = aij;
                         }
                    }
               }
          }

#ifdef DEBUG
          std::cerr << "Khtau=[\n";

          for (octave_idx_type i = 0; i < Khtau.rows(); ++i) {
               for (octave_idx_type j = 0; j < Khtau.columns(); ++j) {
                    std::cerr << std::setw(6) << ElemSpring<ScalarType, MatType>::RealPart(Khtau.xelem(i, j)) << ' ';
               }

               std::cerr << "\n";
          }

          std::cerr << "];\n";



          std::cerr << "Hf=[\n";

          for (octave_idx_type i = 0; i < Hf.rows(); ++i) {
               for (octave_idx_type j = 0; j < Hf.columns(); ++j) {
                    std::cerr << std::setw(6) << Hf.xelem(i, j) << ' ';
               }

               std::cerr << "\n";
          }

          std::cerr << "];\n";


          std::cerr << "Khtau_Hf=[\n";

          for (octave_idx_type i = 0; i < Khtau_Hf.rows(); ++i) {
               for (octave_idx_type j = 0; j < Khtau_Hf.columns(); ++j) {
                    std::cerr << std::setw(6) << ElemSpring<ScalarType, MatType>::RealPart(Khtau_Hf.xelem(i, j)) << ' ';
               }

               std::cerr << "\n";
          }

          std::cerr << "];\n";

          std::cerr << "KHtau=[\n";

          for (octave_idx_type i = 0; i < KHtau.rows(); ++i) {
               for (octave_idx_type j = 0; j < KHtau.columns(); ++j) {
                    std::cerr << std::setw(6) << ElemSpring<ScalarType, MatType>::RealPart(KHtau.xelem(i, j)) << ' ';
               }

               std::cerr << "\n";
          }

          std::cerr << "];\n";
#endif
          return KHtau;
     }
};

template
class ElemContact<double, Matrix>;

template
class ElemContact<std::complex<double>, ComplexMatrix>;

class ElemRBE3: public Element
{
public:
     ElemRBE3(ElementTypes::TypeId eltype, octave_idx_type id, const Matrix& X, const Material* material, const int32NDArray& nodes, const RowVector& omega)
          :Element(eltype, id, X, material, nodes),
           omega(omega) {

          FEM_ASSERT(X.rows() == 6);
          FEM_ASSERT(X.columns() > 1);
          FEM_ASSERT(omega.numel() == nodes.numel() - 1);
     }

     virtual void Assemble(MatrixAss& mat, MeshInfo& info, const DofMap& dof, const FemMatrixType eMatType) const override {
          switch (eMatType) {
          case MAT_STIFFNESS:
          case MAT_STIFFNESS_SYM:
          case MAT_STIFFNESS_SYM_L:
          case MAT_STIFFNESS_FLUID_STRUCT_RE:
               break;

          default:
               return;
          }

          Array<octave_idx_type> ndofidx(dim_vector(nodes.numel() * 6, 1), -1);

          for (octave_idx_type inode = 0; inode < nodes.numel(); ++inode) {
               for (octave_idx_type idof = 0; idof < 6; ++idof) {
                    ndofidx.xelem(inode * 6 + idof) = dof.GetNodeDofIndex(nodes.xelem(inode).value() - 1, DofMap::NDOF_DISPLACEMENT, idof);
               }
          }

          Array<octave_idx_type> edofidx(dim_vector(6, 1));

          for (octave_idx_type idof = 0; idof < edofidx.rows(); ++idof) {
               edofidx.xelem(idof) = dof.GetElemDofIndex(DofMap::ELEM_RBE3, id - 1, idof);
          }

          const octave_idx_type Xrows = X.rows();
          constexpr octave_idx_type xirows = 3;
          const octave_idx_type xicols = nodes.numel() - 1;

          Matrix xi(xirows, xicols);

          for (octave_idx_type j = 1; j < nodes.numel(); ++j) {
               for (octave_idx_type i = 0; i < 3; ++i) {
                    xi.xelem(i + 3 * (j - 1)) = X.xelem(i + Xrows * j) - X.xelem(i);
               }
          }

          FEM_TRACE("xi=[\n" << xi << "];\n");

          const octave_idx_type Srows = xicols * 6;

          Matrix S(Srows, 6);

          for (octave_idx_type k = 0; k < xicols; ++k) {
               for (octave_idx_type j = 0; j < 6; ++j) {
                    for (octave_idx_type i = 0; i < 6; ++i) {
                         const bool alpha = j < 3 || ndofidx.xelem((k + 1) * 6 + j) >= 0;
                         S.xelem(6 * k + i + Srows * j) = alpha * (i == j);
                    }
               }

               /*
                 skew(xi) =
                 0  -z   y
                 z   0  -x
                 -y   x   0

                 -skew(xi) =
                 3   4   5
                 0   z  -y  0
                 -z   0   x  1
                 y  -x   0  2

                 xi =
                 x 0
                 y 1
                 z 2
               */

               S.xelem(6 * k + 0 + Srows * 4) =  xi.xelem(2 + xirows * k);
               S.xelem(6 * k + 0 + Srows * 5) = -xi.xelem(1 + xirows * k);
               S.xelem(6 * k + 1 + Srows * 3) = -xi.xelem(2 + xirows * k);
               S.xelem(6 * k + 1 + Srows * 5) =  xi.xelem(0 + xirows * k);
               S.xelem(6 * k + 2 + Srows * 3) =  xi.xelem(1 + xirows * k);
               S.xelem(6 * k + 2 + Srows * 4) = -xi.xelem(0 + xirows * k);
          }

          FEM_TRACE("S=[\n" << S << "];\n");

          double Lc2 = 0.;

          for (octave_idx_type k = 0; k < xicols; ++k) {
               double norm_xik = 0.;

               for (octave_idx_type i = 0; i < xirows; ++i) {
                    norm_xik += std::pow(xi.xelem(i + xirows * k), 2);
               }

               Lc2 += sqrt(norm_xik);
          }

          Lc2 /= xicols;
          Lc2 *= Lc2;

          FEM_TRACE("Lc2=" << Lc2 << ";\n");

          ColumnVector W(xicols * 6);

          for (octave_idx_type k = 0; k < xicols; ++k) {
               const double omegak = omega.xelem(k);

               for (octave_idx_type i = 0; i < 6; ++i) {
                    W.xelem(k * 6 + i) = omegak * (i < 3 ? 1. : Lc2);
               }
          }

          Matrix STWS(6, 6, 0.);

          for (octave_idx_type j = 0; j < 6; ++j) {
               for (octave_idx_type i = 0; i < 6; ++i) {
                    for (octave_idx_type k = 0; k < S.rows(); ++k) {
                         STWS.xelem(i + 6 * j) += S.xelem(k + Srows * i) * W.xelem(k) * S.xelem(k + Srows * j);
                    }
               }
          }

          FEM_TRACE("f(S^T*W*S)=[\n" << (S.transpose() * DiagMatrix(W) * S - STWS) << "];\n");

          octave_idx_type infoSTWS;

          Matrix X(STWS.inverse(infoSTWS));

          FEM_TRACE("X=[\n" << X << "];\n");

          if (infoSTWS != 0) {
               std::ostringstream os;

               os << "rbe3 element: id " << id << ": X matrix is singular";

               throw std::runtime_error(os.str());
          }

          const octave_idx_type Brows = xicols * 6;

          Matrix B(Brows, 6);

          for (octave_idx_type l = 0; l < xicols; ++l) {
               for (octave_idx_type j = 0; j < 6; ++j) {
                    for (octave_idx_type i = 0; i < 6; ++i) {
                         double Bijl = 0.;

                         for (octave_idx_type k = 0; k < 6; ++k) {
                              Bijl += W.xelem(l * 6 + i) * S.xelem(l * 6 + i + Srows * k) * X.xelem(k + Xrows * j);
                         }

                         B.xelem(l * 6 + i + Brows * j) = Bijl;
                    }
               }
          }

          FEM_TRACE("f=[\n" << (B - DiagMatrix(W) * S * X) << "];\n");
          FEM_TRACE("W=[\n" << W << "];\n");
          FEM_TRACE("B=[\n" << B << "];\n");

          const double beta = mat.GetMatrixInfo(eMatType).beta;

          for (octave_idx_type i = 0; i < 6; ++i) {
               mat.Insert(-beta, ndofidx.xelem(i), edofidx.xelem(i));
               mat.Insert(-beta, edofidx.xelem(i), ndofidx.xelem(i));
          }

          for (octave_idx_type j = 0; j < 6; ++j) {
               for (octave_idx_type i = 0; i < xicols * 6; ++i) {
                    const double Bij = beta * B.xelem(i + Brows * j);
                    mat.Insert(Bij, ndofidx.xelem(i + 6), edofidx.xelem(j));
                    mat.Insert(Bij, edofidx.xelem(j), ndofidx.xelem(i + 6));
               }
          }
     }

     virtual octave_idx_type iGetWorkSpaceSize(FemMatrixType eMatType) const override {
          switch (eMatType) {
          case MAT_STIFFNESS:
          case MAT_STIFFNESS_SYM:
          case MAT_STIFFNESS_SYM_L:
          case MAT_STIFFNESS_FLUID_STRUCT_RE:
               return 8 * 6 + 4 * 6 * 6 * (X.columns() - 1);
          default:
               return 0;
          }
     }

     static constexpr bool bNeedMatrixInfo(Element::FemMatrixType eMatType) {
          switch (eMatType) {
          case MAT_STIFFNESS:
          case MAT_STIFFNESS_FLUID_STRUCT_RE:
          case MAT_STIFFNESS_SYM:
          case MAT_STIFFNESS_SYM_L:
               return true;
          default:
               return false;
          }
     }

     static void AllocIntegrationRule(Element::FemMatrixType eMatType) {
     }
private:
     const RowVector omega;
};

struct BeamCrossSection {
     double A, Ay, Az, It, Iy, Iz;
};

class ElemBeam2: public Element
{
public:
     ElemBeam2(ElementTypes::TypeId eltype, octave_idx_type id, const Matrix& X, const Material* material, const int32NDArray& nodes, const BeamCrossSection& oSect, const ColumnVector& e2, const Matrix& g)
          :Element(eltype, id, X, material, nodes), R(3, 3), g(g) {

          FEM_ASSERT(X.rows() == 6);
          FEM_ASSERT(X.columns() == 2);
          FEM_ASSERT(nodes.numel() == 2);
          FEM_ASSERT(e2.numel() == 3);
          FEM_ASSERT(g.rows() == 3);

          const octave_idx_type Xrows = X.rows();

          l = 0;

          for (octave_idx_type i = 0; i < 3; ++i) {
               double dX = X.xelem(i + Xrows) - X.xelem(i);
               l += dX * dX;
               R.xelem(i) = dX;
          }

          l = sqrt(l);

          if (l == 0) {
               throw std::runtime_error("beam2: zero beam length detected");
          }

          R.xelem(0 + 3 * 2) = R.xelem(1 + 3 * 0) * e2.xelem(2) - R.xelem(2 + 3 * 0) * e2.xelem(1);
          R.xelem(1 + 3 * 2) = R.xelem(2 + 3 * 0) * e2.xelem(0) - R.xelem(0 + 3 * 0) * e2.xelem(2);
          R.xelem(2 + 3 * 2) = R.xelem(0 + 3 * 0) * e2.xelem(1) - R.xelem(1 + 3 * 0) * e2.xelem(0);

          R.xelem(0 + 3 * 1) = R.xelem(1 + 3 * 2) * R.xelem(2 + 3 * 0) - R.xelem(2 + 3 * 2) * R.xelem(1 + 3 * 0);
          R.xelem(1 + 3 * 1) = R.xelem(2 + 3 * 2) * R.xelem(0 + 3 * 0) - R.xelem(0 + 3 * 2) * R.xelem(2 + 3 * 0);
          R.xelem(2 + 3 * 1) = R.xelem(0 + 3 * 2) * R.xelem(1 + 3 * 0) - R.xelem(1 + 3 * 2) * R.xelem(0 + 3 * 0);

          for (octave_idx_type j = 0; j < 3; ++j) {
               double n = 0;

               for (octave_idx_type i = 0; i < 3; ++i) {
                    double Rij = R.xelem(i + 3 * j);
                    n += Rij * Rij;
               }

               n = sqrt(n);

               if (n == 0.) {
                    throw std::runtime_error("beam2: orientation of beam cross section is not valid");
               }

               for (octave_idx_type i = 0; i < 3; ++i) {
                    R.xelem(i + 3 * j) /= n;
               }
          }

          const double E = material->YoungsModulus();
          const double G = material->ShearModulus();
          const double rho = material->Density();

          EA = E * oSect.A;
          GAy = G * oSect.Ay;
          GAz = G * oSect.Az;
          GIt = G * oSect.It;
          EIy = E * oSect.Iy;
          EIz = E * oSect.Iz;
          rhoA = rho * oSect.A;
          rhoIy = rho * oSect.Iy;
          rhoIz = rho * oSect.Iz;
          rhoIp = rhoIy + rhoIz;
          ky = EIz / (l*l * GAy);
          kz = EIy / (l*l * GAz);
     }

     ElemBeam2(const ElemBeam2& oElem)=default;

     virtual ~ElemBeam2()=default;

     virtual void Assemble(MatrixAss& mat, MeshInfo& info, const DofMap& dof, FemMatrixType eMatType) const override {
          void (ElemBeam2::*pfn)(MatrixAss&, const DofMap&) const;

          switch (eMatType) {
          case MAT_STIFFNESS:
          case MAT_STIFFNESS_SYM:
          case MAT_STIFFNESS_SYM_L:
          case MAT_STIFFNESS_FLUID_STRUCT_RE:
               pfn = &ElemBeam2::GlobalStiffnessMatrix;
               break;

          case MAT_STIFFNESS_IM:
          case MAT_STIFFNESS_FLUID_STRUCT_IM:
               pfn = &ElemBeam2::GlobalStiffnessMatrixTanDelta;
               break;

          case MAT_MASS:
          case MAT_MASS_SYM:
          case MAT_MASS_SYM_L:
          case MAT_MASS_FLUID_STRUCT_RE:
               pfn = &ElemBeam2::GlobalMassMatrix;
               break;

          case MAT_DAMPING:
          case MAT_DAMPING_SYM:
          case MAT_DAMPING_SYM_L:
          case MAT_DAMPING_FLUID_STRUCT_RE:
               pfn = &ElemBeam2::GlobalDampingMatrix;
               break;

          case VEC_LOAD_CONSISTENT:
          case VEC_LOAD_LUMPED:
          case VEC_LOAD_FLUID_STRUCT:
               pfn = &ElemBeam2::GlobalLoadVector;
               break;

          default:
               return;
          }

          (this->*pfn)(mat, dof);
     }

     virtual octave_idx_type iGetWorkSpaceSize(FemMatrixType eMatType) const override {
          switch (eMatType) {
          case MAT_STIFFNESS:
          case MAT_STIFFNESS_IM:
          case MAT_STIFFNESS_SYM:
          case MAT_STIFFNESS_SYM_L:
          case MAT_STIFFNESS_FLUID_STRUCT_RE:
          case MAT_STIFFNESS_FLUID_STRUCT_IM:
          case MAT_MASS:
          case MAT_MASS_SYM:
          case MAT_MASS_SYM_L:
          case MAT_MASS_FLUID_STRUCT_RE:
          case MAT_DAMPING:
          case MAT_DAMPING_SYM:
          case MAT_DAMPING_SYM_L:
          case MAT_DAMPING_FLUID_STRUCT_RE:
               return 12 * 12;

          case VEC_LOAD_CONSISTENT:
          case VEC_LOAD_LUMPED:
          case VEC_LOAD_FLUID_STRUCT:
               return 12 * g.columns();

          default:
               return 0;
          }
     }

     virtual double dGetMass() const override {
          return rhoA * l;
     }

     void LocalMatToGlobalMat(Matrix& A) const {
          constexpr octave_idx_type Arows = 12;
          constexpr octave_idx_type Acols = 12;

          FEM_ASSERT(A.rows() == Arows);
          FEM_ASSERT(A.columns() == Acols);
          FEM_ASSERT(R.rows() == 3);
          FEM_ASSERT(R.columns() == 3);

          Matrix RA(Arows, Acols);

          for (octave_idx_type i0 = 0; i0 < 4; ++i0) {
               for (octave_idx_type j0 = 0; j0 < 4; ++j0) {
                    for (octave_idx_type i1 = 0; i1 < 3; ++i1) {
                         for (octave_idx_type j1 = 0; j1 < 3; ++j1) {
                              double RAij = 0;

                              for (octave_idx_type k = 0; k < 3; ++k) {
                                   RAij += R.xelem(i1 + 3 * k) * A.xelem(3 * i0 + k + Arows * (3 * j0 + j1));
                              }

                              RA.xelem(3 * i0 + i1 + Arows * (3 * j0 + j1)) = RAij;
                         }
                    }
               }
          }

          for (octave_idx_type i0 = 0; i0 < 4; ++i0) {
               for (octave_idx_type j0 = 0; j0 < 4; ++j0) {
                    for (octave_idx_type i1 = 0; i1 < 3; ++i1) {
                         for (octave_idx_type j1 = 0; j1 < 3; ++j1) {
                              double Aij = 0;

                              for (octave_idx_type k = 0; k < 3; ++k) {
                                   Aij += R.xelem(j1 + 3 * k) * RA.xelem(3 * i0 + i1 + Arows * (3 * j0 + k));
                              }

                              A.xelem(3 * i0 + i1 + Arows * (3 * j0 + j1)) = Aij;
                         }
                    }
               }
          }

          for (octave_idx_type j = 0; j < Arows; ++j) {
               for (octave_idx_type i = 0; i < j; ++i) {
                    A.xelem(j + Arows * i) = A.xelem(i + Arows * j);
               }
          }
     }

     void LocalVecToGlobalVec(const Matrix& A, Matrix& RA) const {
          constexpr octave_idx_type Arows = 12;
          const octave_idx_type Acols = A.columns();

          FEM_ASSERT(A.rows() == Arows);
          FEM_ASSERT(RA.rows() == Arows);
          FEM_ASSERT(A.columns() == RA.columns());

          for (octave_idx_type j = 0; j < Acols; ++j) {
               for (octave_idx_type i0 = 0; i0 < 4; ++i0) {
                    for (octave_idx_type i1 = 0; i1 < 3; ++i1) {
                         double RAij = 0;

                         for (octave_idx_type k = 0; k < 3; ++k) {
                              RAij += R.xelem(i1 + 3 * k) * A.xelem(3 * i0 + k + Arows * j);
                         }

                         RA.xelem(3 * i0 + i1 + Arows * j) = RAij;
                    }
               }
          }
     }

     void AssembleMatrix(MatrixAss& mat, const Matrix& A, const DofMap& dof) const {
          Array<octave_idx_type> ndofidx(dim_vector(nodes.numel() * 6, 1), -1);

          for (octave_idx_type inode = 0; inode < nodes.numel(); ++inode) {
               for (octave_idx_type idof = 0; idof < 6; ++idof) {
                    ndofidx.xelem(inode * 6 + idof) = dof.GetNodeDofIndex(nodes.xelem(inode).value() - 1, DofMap::NDOF_DISPLACEMENT, idof);
               }
          }

          mat.Insert(A, ndofidx, ndofidx);
     }

     void AssembleVector(MatrixAss& mat, const Matrix& A, const DofMap& dof) const {
          constexpr octave_idx_type nnodes = 2;
#ifdef DEBUG
          constexpr octave_idx_type Arows = nnodes * 6;
#endif
          const octave_idx_type Acols = A.columns();

          FEM_ASSERT(nodes.numel() == nnodes);
          FEM_ASSERT(A.rows() == Arows);

          Array<octave_idx_type> ndofidx(dim_vector(nnodes * 6, 1), -1);

          for (octave_idx_type inode = 0; inode < nnodes; ++inode) {
               for (octave_idx_type idof = 0; idof < 6; ++idof) {
                    ndofidx.xelem(inode * 6 + idof) = dof.GetNodeDofIndex(nodes.xelem(inode).value() - 1, DofMap::NDOF_DISPLACEMENT, idof);
               }
          }

          Array<octave_idx_type> colidx(dim_vector(Acols, 1), -1);

          for (octave_idx_type j = 0; j < Acols; ++j) {
               colidx.xelem(j) = j + 1;
          }

          mat.Insert(A, ndofidx, colidx);
     }

     void GlobalStiffnessMatrix(MatrixAss& mat, const DofMap& dof) const {
          Matrix Ke(12, 12, 0.);

          LocalStiffnessMatrix(Ke);
          LocalMatToGlobalMat(Ke);
          AssembleMatrix(mat, Ke, dof);
     }

     void GlobalStiffnessMatrixTanDelta(MatrixAss& mat, const DofMap& dof) const {
          Matrix Ke_im(12, 12, 0.);

          LocalStiffnessMatrixTanDelta(Ke_im);
          LocalMatToGlobalMat(Ke_im);
          AssembleMatrix(mat, Ke_im, dof);
     }

     void GlobalMassMatrix(MatrixAss& mat, const DofMap& dof) const {
          Matrix Me(12, 12, 0.);

          LocalMassMatrix(Me);
          LocalMatToGlobalMat(Me);
          AssembleMatrix(mat, Me, dof);
     }

     void GlobalDampingMatrix(MatrixAss& mat, const DofMap& dof) const {
          Matrix De(12, 12, 0.);

          LocalDampingMatrix(De);
          LocalMatToGlobalMat(De);
          AssembleMatrix(mat, De, dof);
     }

     void GlobalLoadVector(MatrixAss& mat, const DofMap& dof) const {
          Matrix Re(12, g.columns(), 0.), Rg(12, g.columns());

          LocalLoadVector(Re);
          LocalVecToGlobalVec(Re, Rg);
          AssembleVector(mat, Rg, dof);
     }

     void LocalStiffnessMatrix(Matrix& Ke) const {
          constexpr octave_idx_type Krows = 12;
#ifdef DEBUG
          constexpr octave_idx_type Kcols = 12;
#endif
          FEM_ASSERT(Ke.rows() == Krows);
          FEM_ASSERT(Ke.columns() == Kcols);

          const double l2 = l * l;
          const double l3 = l2 * l;

          Ke.xelem(0 + Krows * 0) = (EA)/l;
          Ke.xelem(1 + Krows * 1) = (12*EIz)/((12*ky+1)*l3);
          Ke.xelem(2 + Krows * 2) = (12*EIy)/((12*kz+1)*l3);
          Ke.xelem(3 + Krows * 3) = (GIt)/l;
          Ke.xelem(2 + Krows * 4) = -(6*EIy)/((12*kz+1)*l2);
          Ke.xelem(4 + Krows * 4) = (4*EIy*(3*kz+1))/((12*kz+1)*l);
          Ke.xelem(1 + Krows * 5) = (6*EIz)/((12*ky+1)*l2);
          Ke.xelem(5 + Krows * 5) = (4*EIz*(3*ky+1))/((12*ky+1)*l);
          Ke.xelem(0 + Krows * 6) = -(EA)/l;
          Ke.xelem(6 + Krows * 6) = (EA)/l;
          Ke.xelem(1 + Krows * 7) = -(12*EIz)/((12*ky+1)*l3);
          Ke.xelem(5 + Krows * 7) = -(6*EIz)/((12*ky+1)*l2);
          Ke.xelem(7 + Krows * 7) = (12*EIz)/((12*ky+1)*l3);
          Ke.xelem(2 + Krows * 8) = -(12*EIy)/((12*kz+1)*l3);
          Ke.xelem(4 + Krows * 8) = (6*EIy)/((12*kz+1)*l2);
          Ke.xelem(8 + Krows * 8) = (12*EIy)/((12*kz+1)*l3);
          Ke.xelem(3 + Krows * 9) = -(GIt)/l;
          Ke.xelem(9 + Krows * 9) = (GIt)/l;
          Ke.xelem(2 + Krows * 10) = -(6*EIy)/((12*kz+1)*l2);
          Ke.xelem(4 + Krows * 10) = -(2*EIy*(6*kz-1))/((12*kz+1)*l);
          Ke.xelem(8 + Krows * 10) = (6*EIy)/((12*kz+1)*l2);
          Ke.xelem(10 + Krows * 10) = (4*EIy*(3*kz+1))/((12*kz+1)*l);
          Ke.xelem(1 + Krows * 11) = (6*EIz)/((12*ky+1)*l2);
          Ke.xelem(5 + Krows * 11) = -(2*EIz*(6*ky-1))/((12*ky+1)*l);
          Ke.xelem(7 + Krows * 11) = -(6*EIz)/((12*ky+1)*l2);
          Ke.xelem(11 + Krows * 11) = (4*EIz*(3*ky+1))/((12*ky+1)*l);
     }

     void LocalMassMatrix(Matrix& Me) const {
          constexpr octave_idx_type Mrows = 12;
#ifdef DEBUG
          constexpr octave_idx_type Mcols = 12;
#endif
          FEM_ASSERT(Me.rows() == Mrows);
          FEM_ASSERT(Me.columns() == Mcols);

          const double l2 = l * l;
          const double ky2 = ky * ky;
          const double kz2 = kz * kz;
          const double b1 = (12*ky+1);
          const double b2 = (12*kz+1);
          const double a1 = b1 * b1;
          const double a2 = b2 * b2;

          Me.xelem(0 + Mrows * 0) = (l*rhoA)/3.0E+0;
          Me.xelem(1 + Mrows * 1) = ((42*rhoIz+1680*ky2*l2*rhoA+294*ky*l2*rhoA+13*l2*rhoA)/(a1*l))/3.5E+1;
          Me.xelem(2 + Mrows * 2) = ((42*rhoIy+1680*kz2*l2*rhoA+294*kz*l2*rhoA+13*l2*rhoA)/(a2*l))/3.5E+1;
          Me.xelem(3 + Mrows * 3) = (l*rhoIp)/3.0E+0;
          Me.xelem(2 + Mrows * 4) = ((1260*kz*rhoIy-21*rhoIy-1260*kz2*l2*rhoA-231*kz*l2*rhoA-11*l2*rhoA)/a2)/2.1E+2;
          Me.xelem(4 + Mrows * 4) = ((l*(5040*kz2*rhoIy+210*kz*rhoIy+14*rhoIy+126*kz2*l2*rhoA+21*kz*l2*rhoA+l2*rhoA))/a2)/1.05E+2;
          Me.xelem(1 + Mrows * 5) = -((1260*ky*rhoIz-21*rhoIz-1260*ky2*l2*rhoA-231*ky*l2*rhoA-11*l2*rhoA)/a1)/2.1E+2;
          Me.xelem(5 + Mrows * 5) = ((l*(5040*ky2*rhoIz+210*ky*rhoIz+14*rhoIz+126*ky2*l2*rhoA+21*ky*l2*rhoA+l2*rhoA))/a1)/1.05E+2;
          Me.xelem(0 + Mrows * 6) = (l*rhoA)/6.0E+0;
          Me.xelem(6 + Mrows * 6) = (l*rhoA)/3.0E+0;
          Me.xelem(1 + Mrows * 7) = ((-3.0E+0)*(28*rhoIz-560*ky2*l2*rhoA-84*ky*l2*rhoA-3*l2*rhoA))/(7.0E+1*a1*l);
          Me.xelem(5 + Mrows * 7) = ((2520*ky*rhoIz-42*rhoIz+2520*ky2*l2*rhoA+378*ky*l2*rhoA+13*l2*rhoA)/a1)/4.2E+2;
          Me.xelem(7 + Mrows * 7) = ((42*rhoIz+1680*ky2*l2*rhoA+294*ky*l2*rhoA+13*l2*rhoA)/(a1*l))/3.5E+1;
          Me.xelem(2 + Mrows * 8) = ((-3.0E+0)*(28*rhoIy-560*kz2*l2*rhoA-84*kz*l2*rhoA-3*l2*rhoA))/(7.0E+1*a2*l);
          Me.xelem(4 + Mrows * 8) = -((2520*kz*rhoIy-42*rhoIy+2520*kz2*l2*rhoA+378*kz*l2*rhoA+13*l2*rhoA)/a2)/4.2E+2;
          Me.xelem(8 + Mrows * 8) = ((42*rhoIy+1680*kz2*l2*rhoA+294*kz*l2*rhoA+13*l2*rhoA)/(a2*l))/3.5E+1;
          Me.xelem(3 + Mrows * 9) = (l*rhoIp)/6.0E+0;
          Me.xelem(9 + Mrows * 9) = (l*rhoIp)/3.0E+0;
          Me.xelem(2 + Mrows * 10) = ((2520*kz*rhoIy-42*rhoIy+2520*kz2*l2*rhoA+378*kz*l2*rhoA+13*l2*rhoA)/a2)/4.2E+2;
          Me.xelem(4 + Mrows * 10) = ((l*(10080*kz2*rhoIy-840*kz*rhoIy-14*rhoIy-504*kz2*l2*rhoA-84*kz*l2*rhoA-3*l2*rhoA))/a2)/4.2E+2;
          Me.xelem(8 + Mrows * 10) = -((1260*kz*rhoIy-21*rhoIy-1260*kz2*l2*rhoA-231*kz*l2*rhoA-11*l2*rhoA)/a2)/2.1E+2;
          Me.xelem(10 + Mrows * 10) = ((l*(5040*kz2*rhoIy+210*kz*rhoIy+14*rhoIy+126*kz2*l2*rhoA+21*kz*l2*rhoA+l2*rhoA))/a2)/1.05E+2;
          Me.xelem(1 + Mrows * 11) = -((2520*ky*rhoIz-42*rhoIz+2520*ky2*l2*rhoA+378*ky*l2*rhoA+13*l2*rhoA)/a1)/4.2E+2;
          Me.xelem(5 + Mrows * 11) = ((l*(10080*ky2*rhoIz-840*ky*rhoIz-14*rhoIz-504*ky2*l2*rhoA-84*ky*l2*rhoA-3*l2*rhoA))/a1)/4.2E+2;
          Me.xelem(7 + Mrows * 11) = ((1260*ky*rhoIz-21*rhoIz-1260*ky2*l2*rhoA-231*ky*l2*rhoA-11*l2*rhoA)/a1)/2.1E+2;
          Me.xelem(11 + Mrows * 11) = ((l*(5040*ky2*rhoIz+210*ky*rhoIz+14*rhoIz+126*ky2*l2*rhoA+21*ky*l2*rhoA+l2*rhoA))/a1)/1.05E+2;
     }

     void LocalDampingMatrix(Matrix& De) const {
          const double alpha = material->AlphaDamping();
          const double beta = material->BetaDamping();

          constexpr octave_idx_type Drows = 12;
          constexpr octave_idx_type Dcols = 12;

          Matrix Me(Drows, Dcols, 0.);

          LocalMassMatrix(Me);
          LocalStiffnessMatrix(De);

          for (octave_idx_type j = 0; j < Dcols; ++j) {
               for (octave_idx_type i = 0; i < Drows; ++i) {
                    De.xelem(i + Drows * j) = alpha * Me.xelem(i + Drows * j) + beta * De.xelem(i + Drows * j);
               }
          }
     }

     void LocalStiffnessMatrixTanDelta(Matrix& Ke_im) const {
          const double tan_delta = material->TanDeltaDamping();
          constexpr octave_idx_type Krows = 12, Kcols = 12;

          LocalStiffnessMatrix(Ke_im);

          for (octave_idx_type j = 0; j < Kcols; ++j) {
               for (octave_idx_type i = 0; i < Krows; ++i) {
                    Ke_im.xelem(i + Krows * j) *= tan_delta;
               }
          }
     }

     void LocalLoadVector(Matrix& Pe) const {
          constexpr octave_idx_type Rrows = 3;
          constexpr octave_idx_type Prows = 12;
          constexpr octave_idx_type grows = 3;
          const octave_idx_type gcols = g.columns();

          FEM_ASSERT(R.rows() == Rrows);
          FEM_ASSERT(R.columns() == Rrows);
          FEM_ASSERT(g.rows() == grows);
          FEM_ASSERT(Pe.columns() == g.columns());
          FEM_ASSERT(Pe.rows() == Prows);

          const double l2 = l * l;
          ColumnVector pe(grows);

          for (octave_idx_type i = 0; i < gcols; ++i) {
               pe.xelem(0) = (g.xelem(2 + grows * i) * R.xelem(2 + Rrows * 0) + g.xelem(1 + grows * i) * R.xelem(1 + Rrows * 0) + g.xelem(0 + grows * i) * R.xelem(0 + Rrows * 0)) * rhoA;
               pe.xelem(1) = (g.xelem(2 + grows * i) * R.xelem(2 + Rrows * 1) + g.xelem(1 + grows * i) * R.xelem(1 + Rrows * 1) + g.xelem(0 + grows * i) * R.xelem(0 + Rrows * 1)) * rhoA;
               pe.xelem(2) = (g.xelem(2 + grows * i) * R.xelem(2 + Rrows * 2) + g.xelem(1 + grows * i) * R.xelem(1 + Rrows * 2) + g.xelem(0 + grows * i) * R.xelem(0 + Rrows * 2)) * rhoA;

               Pe.xelem(0 + Prows * i) = (pe.xelem(0) * l) / 2.0E+0;
               Pe.xelem(1 + Prows * i) = (pe.xelem(1) * l) / 2.0E+0;
               Pe.xelem(2 + Prows * i) = (pe.xelem(2) * l) / 2.0E+0;
               Pe.xelem(3 + Prows * i) = 0;
               Pe.xelem(4 + Prows * i) = -(pe.xelem(2) * l2) / 1.2E+1;
               Pe.xelem(5 + Prows * i) = (pe.xelem(1) * l2) / 1.2E+1;
               Pe.xelem(6 + Prows * i) = (pe.xelem(0) * l) / 2.0E+0;
               Pe.xelem(7 + Prows * i) = (pe.xelem(1) * l) / 2.0E+0;
               Pe.xelem(8 + Prows * i) = (pe.xelem(2) * l) / 2.0E+0;
               Pe.xelem(9 + Prows * i) = 0;
               Pe.xelem(10 + Prows * i) = (pe.xelem(2) * l2) / 1.2E+1;
               Pe.xelem(11 + Prows * i) = -(pe.xelem(1) * l2) / 1.2E+1;
          }
     }

     static void AllocIntegrationRule(Element::FemMatrixType eMatType) {
     }
private:
     Matrix R;
     const Matrix g;
     double l, ky, kz, EA, GAy, GAz, GIt, EIy, EIz, rhoA, rhoIp, rhoIy, rhoIz;
};

class ElemRigidBody: public Element
{
public:
     ElemRigidBody(ElementTypes::TypeId eltype, octave_idx_type id, const Matrix& X, const Material* material, const int32NDArray& nodes, double m, const Matrix& J, const ColumnVector& lcg, const Matrix& g)
          :Element(eltype, id, X, material, nodes), m(m), J(J), skew_lcg(3, 3, 0.), g(g) {

          skew_lcg.xelem(0 + 3 * 1) = -lcg.xelem(2);
          skew_lcg.xelem(1 + 3 * 0) = lcg.xelem(2);
          skew_lcg.xelem(0 + 3 * 2) = lcg.xelem(1);
          skew_lcg.xelem(2 + 3 * 0) = -lcg.xelem(1);
          skew_lcg.xelem(1 + 3 * 2) = -lcg.xelem(0);
          skew_lcg.xelem(2 + 3 * 1) = lcg.xelem(0);

          FEM_ASSERT(X.rows() == 6);
          FEM_ASSERT(X.columns() == 1);
          FEM_ASSERT(nodes.numel() == 1);
     }

     ElemRigidBody(const ElemRigidBody& oElem)=default;

     virtual ~ElemRigidBody()=default;

     virtual void Assemble(MatrixAss& mat, MeshInfo& info, const DofMap& dof, FemMatrixType eMatType) const override {
          void (ElemRigidBody::*pfn)(MatrixAss& mat, const DofMap& dof) const;

          switch (eMatType) {
          case MAT_MASS:
          case MAT_MASS_SYM:
          case MAT_MASS_SYM_L:
          case MAT_MASS_FLUID_STRUCT_RE:
               pfn = &ElemRigidBody::GlobalMassMatrix;
               break;

          case VEC_LOAD_CONSISTENT:
          case VEC_LOAD_LUMPED:
          case VEC_LOAD_FLUID_STRUCT:
               pfn = &ElemRigidBody::GlobalLoadVector;
               break;

          default:
               return;
          }

          (this->*pfn)(mat, dof);
     }

     virtual octave_idx_type iGetWorkSpaceSize(FemMatrixType eMatType) const override {
          switch (eMatType) {
          case MAT_MASS:
          case MAT_MASS_SYM:
          case MAT_MASS_SYM_L:
          case MAT_MASS_FLUID_STRUCT_RE:
               return 6 * 6;

          case VEC_LOAD_CONSISTENT:
          case VEC_LOAD_LUMPED:
          case VEC_LOAD_FLUID_STRUCT:
               return 6 * g.columns();

          default:
               return 0;
          }
     }

     virtual double dGetMass() const override {
          return m;
     }

     void GlobalMassMatrix(MatrixAss& mat, const DofMap& dof) const {
          Matrix A(6, 6, 0.);

          LocalMassMatrix(A);

          Array<octave_idx_type> ndofidx(dim_vector(nodes.numel() * 6, 1), -1);

          const octave_idx_type inode = nodes.xelem(0).value() - 1;

          for (octave_idx_type idof = 0; idof < 6; ++idof) {
               ndofidx.xelem(idof) = dof.GetNodeDofIndex(inode, DofMap::NDOF_DISPLACEMENT, idof);
          }

          mat.Insert(A, ndofidx, ndofidx);
     }

     void GlobalLoadVector(MatrixAss& mat, const DofMap& dof) const {
          Matrix Re(6, g.columns());

          LocalLoadVector(Re);

          Array<octave_idx_type> ndofidx(dim_vector(nodes.numel() * 6, 1), -1);

          const octave_idx_type inode = nodes.xelem(0).value() - 1;

          for (octave_idx_type idof = 0; idof < 6; ++idof) {
               ndofidx.xelem(idof) = dof.GetNodeDofIndex(inode, DofMap::NDOF_DISPLACEMENT, idof);
          }

          Array<octave_idx_type> colidx(dim_vector(g.columns(), 1), -1);

          for (octave_idx_type j = 0; j < g.columns(); ++j) {
               colidx.xelem(j) = j + 1;
          }

          mat.Insert(Re, ndofidx, colidx);
     }

     void LocalMassMatrix(Matrix& Me) const {
          constexpr octave_idx_type Mrows = 6;
#ifdef DEBUG
          constexpr octave_idx_type Mcols = 6;
#endif
          FEM_ASSERT(Me.rows() == Mrows);
          FEM_ASSERT(Me.columns() == Mcols);

          for (octave_idx_type i = 0; i < 3; ++i) {
               Me.xelem(i + Mrows * i) = m;
          }

          for (octave_idx_type j = 0; j < 3; ++j) {
               for (octave_idx_type i = 0; i < 3; ++i) {
                    double Jij = J.xelem(i + 3 * j);

                    for (octave_idx_type k = 0; k < 3; ++k) {
                         Jij += m * skew_lcg.xelem(k + 3 * i) * skew_lcg.xelem(k + 3 * j);
                    }

                    Me.xelem(i + 3 + Mrows * (j + 3)) = Jij;
               }
          }

          for (octave_idx_type j = 0; j < 3; ++j) {
               for (octave_idx_type i = 0; i < 3; ++i) {
                    Me.xelem(j + 3 + Mrows * i) = Me.xelem(i + Mrows * (j + 3)) = -m * skew_lcg.xelem(i + 3 * j);
               }
          }
     }

     void LocalLoadVector(Matrix& Re) const {
          FEM_ASSERT(Re.rows() == 6);
          FEM_ASSERT(Re.columns() == g.columns());

          for (octave_idx_type j = 0; j < g.columns(); ++j) {
               for (octave_idx_type i = 0; i < 3; ++i) {
                    Re.xelem(i + 6 * j) = g.xelem(i + 3 * j) * m;
               }

               for (octave_idx_type i = 0; i < 3; ++i) {
                    double Mij = 0.;

                    for (octave_idx_type k = 0; k < 3; ++k) {
                         Mij += skew_lcg.xelem(i + 3 * k) * Re.xelem(k + 6 * j);
                    }

                    Re.xelem(i + 3 + 6 * j) = Mij;
               }
          }
     }

     static void AllocIntegrationRule(Element::FemMatrixType) {
     }
private:
     double m;
     Matrix J;
     Matrix skew_lcg;
     const Matrix g;
};

class Element3D: public Element
{
public:
     struct ElementData {
          ElementData(const octave_map& load_case, const Matrix& nodes, const octave_scalar_map& elements)
               :oRefStrain{load_case, nodes},
                oPreStress{load_case, nodes},
                oGravity{load_case},
                oCentrifugal{load_case},
                oPML{elements} {
          }

          const StrainField oRefStrain;
          const PreStressField oPreStress;
          const GravityLoad oGravity;
          const CentrifugalLoad oCentrifugal;
          PerfectlyMatchedLayer oPML;
     };

     Element3D(ElementTypes::TypeId eltype, octave_idx_type id, const Matrix& X, const Material* material, const int32NDArray& nodes, const ElementData& data)
          :Element(eltype, id, X, material, nodes), eMaterial(material->GetMaterialType()), oElemData(data) {

          FEM_ASSERT(X.rows() == 3);

          const octave_idx_type iNumNodes = nodes.numel();
          const octave_idx_type iNumLoadsTemp = data.oRefStrain.rgTemperature.numel();
          const octave_idx_type iNumLoadsStrain = data.oRefStrain.rgRefStrain.numel();
          const octave_idx_type iNumPreStress = data.oPreStress.rgPreStress.numel();

          iNumPreLoads = std::max(std::max(std::max(iNumLoadsTemp, iNumLoadsStrain), data.oGravity.g.columns()), data.oCentrifugal.WxWx.ndims() >= 3 ? data.oCentrifugal.WxWx.dim3() : 1);

          FEM_ASSERT(iNumLoadsTemp && iNumLoadsStrain ? iNumLoadsTemp == iNumLoadsStrain : true);

          if (iNumLoadsTemp) {
               dTheta.resize(iNumNodes, iNumLoadsTemp);

               for (octave_idx_type j = 0; j < iNumLoadsTemp; ++j) {
                    const NDArray dThetaj = data.oRefStrain.rgTemperature.xelem(j).array_value();

                    for (octave_idx_type i = 0; i < iNumNodes; ++i) {
                         dTheta.xelem(i + iNumNodes * j) = dThetaj.xelem(nodes.xelem(i).value() - 1);
                    }
               }
          }

          if (iNumLoadsStrain) {
               const octave_idx_type iNumStrains = material->LinearElasticity().rows();

               epsilonRef.resize(dim_vector(iNumStrains, iNumNodes, iNumLoadsStrain), 0.);

               for (octave_idx_type k = 0; k < iNumLoadsStrain; ++k) {
                    const octave_scalar_map maEpsilonRefk = data.oRefStrain.rgRefStrain.xelem(k).scalar_map_value();

                    const std::string strElemName = ElementTypes::GetType(eltype).name;
                    const auto iterEpsilonRefk = maEpsilonRefk.seek(strElemName);

                    if (iterEpsilonRefk == maEpsilonRefk.end()) {
                         continue;
                    }

                    const NDArray epsilonRefk = maEpsilonRefk.contents(iterEpsilonRefk).array_value();
                    const octave_idx_type iNumElem = epsilonRefk.dim1();

                    if (epsilonRefk.ndims() != 3 || iNumElem < id || epsilonRefk.dim2() != iNumNodes || ((epsilonRefk.ndims() >= 3 ? epsilonRefk.dim3() : 1) != iNumStrains)) {
                         throw std::runtime_error("fem_ass_matrix: invalid number of dimensions for load_case.epsilon0."s + strElemName);
                    }

                    for (octave_idx_type j = 0; j < iNumNodes; ++j) {
                         for (octave_idx_type i = 0; i < iNumStrains; ++i) {
                              epsilonRef.xelem(i + iNumStrains * (j + iNumNodes * k)) = epsilonRefk.xelem(id - 1 + iNumElem * (j + iNumNodes * i));
                         }
                    }
               }
          }

          if (iNumPreStress) {
               const octave_idx_type iNumStressComp = material->LinearElasticity().rows();

               tauRef.resize(dim_vector(iNumStressComp, iNumNodes, iNumPreStress), 0.);

               for (octave_idx_type k = 0; k < iNumPreStress; ++k) {
                    const octave_scalar_map maTauRefk = data.oPreStress.rgPreStress.xelem(k).scalar_map_value();

                    const std::string strElemName = ElementTypes::GetType(eltype).name;
                    const auto iterTauRefk = maTauRefk.seek(strElemName);

                    if (iterTauRefk == maTauRefk.end()) {
                         continue;
                    }

                    const NDArray tauRefk = maTauRefk.contents(iterTauRefk).array_value();
                    const octave_idx_type iNumElem = tauRefk.dim1();

                    if (tauRefk.ndims() != 3 || iNumElem < id || tauRefk.dim2() != iNumNodes || ((tauRefk.ndims() >= 3 ? tauRefk.dim3() : 1) != iNumStressComp)) {
                         throw std::runtime_error("fem_ass_matrix: invalid number of dimensions for load_case.tau0."s + strElemName);
                    }

                    for (octave_idx_type j = 0; j < iNumNodes; ++j) {
                         for (octave_idx_type i = 0; i < iNumStressComp; ++i) {
                              tauRef.xelem(i + iNumStressComp * (j + iNumNodes * k)) = tauRefk.xelem(id - 1 + iNumElem * (j + iNumNodes * i));
                         }
                    }
               }
          }

          if (!data.oPML.rgElem[eltype].f.isempty()) {
               FEM_ASSERT(data.oPML.rgElem[eltype].f.ndims() == 3);
               FEM_ASSERT(data.oPML.rgElem[eltype].e1.ndims() == 3);
               FEM_ASSERT(data.oPML.rgElem[eltype].e2.ndims() == 3);
               FEM_ASSERT(data.oPML.rgElem[eltype].f.rows() == 3);
               FEM_ASSERT(data.oPML.rgElem[eltype].f.columns() == iNumNodes);
               FEM_ASSERT(data.oPML.rgElem[eltype].f.pages() >= id);
               FEM_ASSERT(data.oPML.rgElem[eltype].e1.rows() == 3);
               FEM_ASSERT(data.oPML.rgElem[eltype].e1.columns() == iNumNodes);
               FEM_ASSERT(data.oPML.rgElem[eltype].e2.rows() == 3);
               FEM_ASSERT(data.oPML.rgElem[eltype].e2.columns() == iNumNodes);
               FEM_ASSERT(data.oPML.rgElem[eltype].e1.pages() >= id);
               FEM_ASSERT(data.oPML.rgElem[eltype].e2.pages() >= id);

               f = data.oPML.rgElem[eltype].f.linear_slice(3 * iNumNodes * (id - 1), 3 * iNumNodes * id).reshape(3, iNumNodes);
               e1 = data.oPML.rgElem[eltype].e1.linear_slice(3 * iNumNodes * (id - 1), 3 * iNumNodes * id).reshape(3, iNumNodes);
               e2 = data.oPML.rgElem[eltype].e2.linear_slice(3 * iNumNodes * (id - 1), 3 * iNumNodes * id).reshape(3, iNumNodes);
          }
     }

     virtual void Assemble(MatrixAss& mat, MeshInfo& info, const DofMap& dof, FemMatrixType eMatType) const override final {
          void (Element3D::*pFunc)(Matrix&, MeshInfo&, FemMatrixType) const;

          const octave_idx_type iNumDof = iGetNumDof(eMatType);

          octave_idx_type iNumRows = 0, iNumCols = 0;

          DofMap::NodalDofType eDofType;

          switch (eMaterial) {
          case Material::MAT_TYPE_SOLID:
               eDofType = DofMap::NDOF_DISPLACEMENT;
               break;
          case Material::MAT_TYPE_THERMAL:
               eDofType = DofMap::NDOF_TEMPERATURE;
               break;
          case Material::MAT_TYPE_FLUID:
               eDofType = DofMap::NDOF_VELOCITY_POT;
               break;
          default:
               throw std::logic_error("element 3D: unknown material type");
          }

          switch (eMatType) {
          case MAT_STIFFNESS:
          case MAT_STIFFNESS_SYM:
          case MAT_STIFFNESS_SYM_L:
               pFunc = &Element3D::StiffnessMatrix;
               iNumRows = iNumCols = iNumDof;
               break;

          case MAT_STIFFNESS_IM:
               pFunc = &Element3D::StiffnessMatrixTanDelta;
               iNumRows = iNumCols = iNumDof;
               break;

          case MAT_STIFFNESS_TAU0:
               pFunc = &Element3D::StiffnessMatrixNL;
               iNumRows = iNumCols = iNumDof;
               break;

          case MAT_STIFFNESS_OMEGA:
               if (oElemData.oCentrifugal.WxWx.numel()) {
                    pFunc = &Element3D::CentrifugalStiffnessMatrix;
                    iNumRows = iNumCols = iNumDof;
               }
               break;

          case MAT_STIFFNESS_OMEGA_DOT:
               if (oElemData.oCentrifugal.WPx.numel()) {
                    pFunc = &Element3D::AngularAccelerationStiffnessMatrix;
                    iNumRows = iNumCols = iNumDof;
               }
               break;

          case MAT_DAMPING_OMEGA:
               if (oElemData.oCentrifugal.Wx.numel()) {
                    pFunc = &Element3D::CoriolisDampingMatrix;
                    iNumRows = iNumCols = iNumDof;
               }
               break;

          case MAT_MASS:
          case MAT_MASS_SYM:
          case MAT_MASS_SYM_L:
          case MAT_MASS_LUMPED:
               pFunc = &Element3D::MassMatrix;
               iNumRows = iNumCols = iNumDof;
               break;

          case MAT_DAMPING:
          case MAT_DAMPING_SYM:
          case MAT_DAMPING_SYM_L:
               pFunc = &Element3D::DampingMatrix;
               iNumRows = iNumCols = iNumDof;
               break;

          case VEC_LOAD_CONSISTENT:
          case VEC_LOAD_LUMPED:
          case VEC_LOAD_FLUID_STRUCT:
               if (eMaterial == Material::MAT_TYPE_SOLID) {
                    pFunc = &Element3D::StructuralLoadVector;
                    iNumRows = iNumDof;
                    iNumCols = iNumPreLoads;
               }
               break;

          case MAT_THERMAL_COND:
               pFunc = &Element3D::ThermalConductivityMatrix;
               iNumRows = iNumCols = iNumDof;
               break;

          case MAT_HEAT_CAPACITY:
               pFunc = &Element3D::HeatCapacityMatrix;
               iNumRows = iNumCols = iNumDof;
               break;

          case MAT_STIFFNESS_ACOUSTICS_RE:
               if (f.isempty()) {
                    pFunc = &Element3D::AcousticStiffnessMatrix;
               } else {
                    pFunc = &Element3D::AcousticStiffnessMatrixPML;
               }

               iNumRows = iNumCols = iNumDof;
               break;

          case MAT_STIFFNESS_ACOUSTICS_IM:
               if (!f.isempty()) {
                    pFunc = &Element3D::AcousticStiffnessMatrixPML;
                    iNumRows = iNumCols = iNumDof;
               }
               break;

          case MAT_DAMPING_ACOUSTICS_RE:
               if (f.isempty()) {
                    pFunc = &Element3D::AcousticDampingMatrix;
               } else {
                    pFunc = &Element3D::AcousticDampingMatrixPML;
               }

               iNumRows = iNumCols = iNumDof;
               break;

          case MAT_DAMPING_ACOUSTICS_IM:
               if (!f.isempty()) {
                    pFunc = &Element3D::AcousticDampingMatrixPML;
                    iNumRows = iNumCols = iNumDof;
               }
               break;

          case MAT_MASS_ACOUSTICS_RE:
               if (f.isempty()) {
                    pFunc = &Element3D::AcousticMassMatrix;
               } else {
                    pFunc = &Element3D::AcousticMassMatrixPML;
               }
               iNumRows = iNumCols = iNumDof;
               break;

          case MAT_MASS_ACOUSTICS_IM:
               if (!f.isempty()) {
                    pFunc = &Element3D::AcousticMassMatrixPML;
                    iNumRows = iNumCols = iNumDof;
               }
               break;

          case MAT_STIFFNESS_FLUID_STRUCT_RE:
               switch (eMaterial) {
               case Material::MAT_TYPE_SOLID:
                    pFunc = &Element3D::StiffnessMatrix;
                    break;
               case Material::MAT_TYPE_FLUID:
                    if (f.isempty()) {
                         pFunc = &Element3D::AcousticStiffnessMatrix;
                    } else {
                         pFunc = &Element3D::AcousticStiffnessMatrixPML;
                    }
                    break;
               default:
                    throw std::logic_error("element 3D: material not supported");
               }

               iNumRows = iNumCols = iNumDof;
               break;

          case MAT_STIFFNESS_FLUID_STRUCT_IM:
               switch (eMaterial) {
               case Material::MAT_TYPE_FLUID:
                    if (!f.isempty()) {
                         pFunc = &Element3D::AcousticStiffnessMatrixPML;
                         iNumRows = iNumCols = iNumDof;
                    }
                    break;
               case Material::MAT_TYPE_SOLID:
                    pFunc = &Element3D::StiffnessMatrixTanDelta;
                    iNumRows = iNumCols = iNumDof;
                    break;
               default:
                    throw std::logic_error("element 3D: material not supported");
               }
               break;

          case MAT_DAMPING_FLUID_STRUCT_RE:
               switch (eMaterial) {
               case Material::MAT_TYPE_SOLID:
                    pFunc = &Element3D::DampingMatrix;
                    break;
               case Material::MAT_TYPE_FLUID:
                    if (f.isempty()) {
                         pFunc = &Element3D::AcousticDampingMatrix;
                    } else {
                         pFunc = &Element3D::AcousticDampingMatrixPML;
                    }
                    break;
               default:
                    throw std::logic_error("element 3D: material not supported");
               }

               iNumRows = iNumCols = iNumDof;
               break;

          case MAT_DAMPING_FLUID_STRUCT_IM:
               if (eMaterial == Material::MAT_TYPE_FLUID && !f.isempty()) {
                    pFunc = &Element3D::AcousticDampingMatrixPML;
                    iNumRows = iNumCols = iNumDof;
               }
               break;

          case MAT_MASS_FLUID_STRUCT_RE:
               switch (eMaterial) {
               case Material::MAT_TYPE_SOLID:
                    pFunc = &Element3D::MassMatrix;
                    break;
               case Material::MAT_TYPE_FLUID:
                    if (f.isempty()) {
                         pFunc = &Element3D::AcousticMassMatrix;
                    } else {
                         pFunc = &Element3D::AcousticMassMatrixPML;
                    }
                    break;
               default:
                    throw std::logic_error("element 3D: material not supported");
               }

               iNumRows = iNumCols = iNumDof;
               break;
          case MAT_MASS_FLUID_STRUCT_IM:
               if (eMaterial == Material::MAT_TYPE_FLUID && !f.isempty()) {
                    pFunc = &Element3D::AcousticMassMatrixPML;
                    iNumRows = iNumCols = iNumDof;
               }
               break;

          default:
               break;
          }

          if (iNumCols == 0) {
               return;
          }

          Array<octave_idx_type> dofidx(dim_vector(iNumDof, 1), 0);

          constexpr unsigned uScalarFieldMask = DofMap::DO_THERMAL | DofMap::DO_ACOUSTICS;
          const octave_idx_type inodemaxdof = ((eMatType & uScalarFieldMask) != 0u || eMaterial == Material::MAT_TYPE_FLUID) ? 1 : 3;
          const octave_idx_type nnodes = nodes.numel();

          for (octave_idx_type inode = 0; inode < nnodes; ++inode) {
               const octave_idx_type inodeidx = nodes.xelem(inode).value() - 1;
               for (octave_idx_type idof = 0; idof < inodemaxdof; ++idof) {
                    dofidx.xelem(inode * inodemaxdof + idof) = dof.GetNodeDofIndex(inodeidx, eDofType, idof);
               }
          }

          Matrix Ae(iNumRows, iNumCols, 0.);

          (this->*pFunc)(Ae, info, eMatType);

          switch (eMatType) {
          case VEC_LOAD_CONSISTENT:
          case VEC_LOAD_LUMPED:
          case VEC_LOAD_FLUID_STRUCT: {
               Array<octave_idx_type> dofidxcol(dim_vector(iNumCols, 1));

               for (octave_idx_type i = 0; i < iNumCols; ++i) {
                    dofidxcol.xelem(i) = i + 1;
               }

               mat.Insert(Ae, dofidx, dofidxcol);
          } break;
          default:
               mat.Insert(Ae, dofidx, dofidx);
          }
     }

     octave_idx_type iGetNumDof(FemMatrixType eMatType) const {
          switch (material->GetMaterialType()) {
          case Material::MAT_TYPE_SOLID:
               return nodes.numel() * 3;
          case Material::MAT_TYPE_THERMAL:
          case Material::MAT_TYPE_FLUID:
               return nodes.numel();
          default:
               throw std::logic_error("element 3D: material not supported");
          }
     }

     virtual octave_idx_type iGetWorkSpaceSize(FemMatrixType eMatType) const override final {
          switch (eMatType) {
          case MAT_MASS:
          case MAT_MASS_FLUID_STRUCT_RE:
          case MAT_MASS_FLUID_STRUCT_IM:
          case MAT_STIFFNESS:
          case MAT_STIFFNESS_IM:
          case MAT_STIFFNESS_TAU0:
          case MAT_STIFFNESS_OMEGA:
          case MAT_STIFFNESS_OMEGA_DOT:
          case MAT_DAMPING_OMEGA:
          case MAT_STIFFNESS_FLUID_STRUCT_RE:
          case MAT_STIFFNESS_FLUID_STRUCT_IM:
          case MAT_DAMPING:
          case MAT_DAMPING_FLUID_STRUCT_RE:
          case MAT_THERMAL_COND:
          case MAT_HEAT_CAPACITY:
          case MAT_MASS_ACOUSTICS_RE:
          case MAT_MASS_ACOUSTICS_IM:
          case MAT_STIFFNESS_ACOUSTICS_RE:
          case MAT_STIFFNESS_ACOUSTICS_IM:
          case MAT_DAMPING_ACOUSTICS_RE:
          case MAT_MASS_SYM:
          case MAT_MASS_SYM_L:
          case MAT_STIFFNESS_SYM:
          case MAT_STIFFNESS_SYM_L:
          case MAT_DAMPING_SYM:
          case MAT_DAMPING_SYM_L:
          case MAT_MASS_LUMPED: {
               const octave_idx_type iNumDof = iGetNumDof(eMatType);

               return iNumDof * iNumDof;
          }
          case VEC_LOAD_CONSISTENT:
          case VEC_LOAD_LUMPED:
          case VEC_LOAD_FLUID_STRUCT:
               return iGetNumDof(eMatType) * iNumPreLoads;

          default:
               return 0;
          }
     }

     virtual void PostProcElem(FemMatrixType eMatType, PostProcData& oSolution) const override final {
          switch (eMaterial) {
          case Material::MAT_TYPE_SOLID:
               switch (eMatType) {
               case VEC_INERTIA_M1:
                    InertiaMoment1(oSolution.GetField(PostProcData::VEC_G_STRUCT_INERTIA_M1_RE, eltype), eMatType);
                    break;

               case MAT_INERTIA_J:
                    InertiaMatrix(oSolution.GetField(PostProcData::MAT_G_STRUCT_INERTIA_J_RE, eltype), eMatType);
                    break;

               case MAT_INERTIA_INV3:
                    InertiaInv3(oSolution.GetField(PostProcData::MAT_G_STRUCT_INERTIA_INV3_RE, eltype), eMatType, oSolution.GetField(PostProcData::VEC_NO_STRUCT_DISPLACEMENT_RE, eltype));
                    break;

               case MAT_INERTIA_INV4:
                    InertiaInv4(oSolution.GetField(PostProcData::MAT_G_STRUCT_INERTIA_INV4_RE, eltype), eMatType, oSolution.GetField(PostProcData::VEC_NO_STRUCT_DISPLACEMENT_RE, eltype));
                    break;

               case MAT_INERTIA_INV5:
                    InertiaInv5(oSolution.GetField(PostProcData::MAT_G_STRUCT_INERTIA_INV5_RE, eltype), eMatType, oSolution.GetField(PostProcData::VEC_NO_STRUCT_DISPLACEMENT_RE, eltype));
                    break;

               case MAT_INERTIA_INV8:
                    InertiaInv8(oSolution.GetField(PostProcData::MAT_G_STRUCT_INERTIA_INV8_RE, eltype), eMatType, oSolution.GetField(PostProcData::VEC_NO_STRUCT_DISPLACEMENT_RE, eltype));
                    break;

               case MAT_INERTIA_INV9:
                    InertiaInv9(oSolution.GetField(PostProcData::MAT_G_STRUCT_INERTIA_INV9_RE, eltype), eMatType, oSolution.GetField(PostProcData::VEC_NO_STRUCT_DISPLACEMENT_RE, eltype));
                    break;

               case VEC_STRESS_CAUCH:
                    StressNodalElem(oSolution.GetField(PostProcData::VEC_EL_STRUCT_STRESS_CAUCH_RE, eltype), eMatType, oSolution.GetField(PostProcData::VEC_NO_STRUCT_DISPLACEMENT_RE, eltype));
                    break;

               case VEC_STRAIN_TOTAL:
                    StrainNodalElem(oSolution.GetField(PostProcData::VEC_EL_STRUCT_STRAIN_TOTAL_RE, eltype), eMatType, oSolution.GetField(PostProcData::VEC_NO_STRUCT_DISPLACEMENT_RE, eltype));
                    break;
               default:
                    break;
               }
               break;
          case Material::MAT_TYPE_FLUID:
               switch (eMatType) {
               case VEC_PARTICLE_VELOCITY:
                    ParticleVelocityNodalElem<double>(oSolution.GetField(PostProcData::VEC_EL_ACOUSTIC_PART_VEL_RE, eltype),
                                                      eMatType,
                                                      oSolution.GetField(PostProcData::SCA_NO_ACOUSTIC_PART_VEL_POT_RE, eltype),
                                                      oSolution.GetField(PostProcData::SCA_NO_ACOUSTIC_PART_VEL_POT_P_RE, eltype));
                    break;

               case VEC_PARTICLE_VELOCITY_C:
                    ParticleVelocityNodalElem<std::complex<double> >(oSolution.GetField(PostProcData::VEC_EL_ACOUSTIC_PART_VEL_C, eltype),
                                                                     eMatType,
                                                                     oSolution.GetField(PostProcData::SCA_NO_ACOUSTIC_PART_VEL_POT_C, eltype),
                                                                     oSolution.GetField(PostProcData::SCA_NO_ACOUSTIC_PART_VEL_POT_P_C, eltype));
                    break;
               default:
                    break;
               }
               break;
          case Material::MAT_TYPE_THERMAL:
               break;
          }

          if (eMatType & MAT_COLL_PNT_OUTPUT) {
               auto eMatTypeColloc = static_cast<FemMatrixType>(eMatType & ~MAT_COLL_PNT_OUTPUT);
               CollocationPoints(oSolution.GetField(PostProcData::VEC_EL_COLLOC_POINTS_RE, eltype), eMatTypeColloc);
          }
     }

     double dGetVolume() const {
          const IntegrationRule& oIntegRule = GetIntegrationRule(MAT_MASS);
          const octave_idx_type iNumDir = oIntegRule.iGetNumDirections();
          ColumnVector rv(iNumDir);
          Matrix J(iNumDir, iNumDir);

          double dV = 0.;

          for (octave_idx_type i = 0; i < oIntegRule.iGetNumEvalPoints(); ++i) {
               const double alpha = oIntegRule.dGetWeight(i);

               for (octave_idx_type j = 0; j < iNumDir; ++j) {
                    rv.xelem(j) = oIntegRule.dGetPosition(i, j);
               }

               const double detJ = Jacobian(rv, J);

               dV += alpha * detJ;
          }

          return dV;
     }

     virtual double dGetMass() const override final {
          return dGetVolume() * material->Density();
     }

     virtual octave_idx_type iGetNumCollocPoints(FemMatrixType eMatType) const override {
          const auto eMatTypeColloc = static_cast<FemMatrixType>(eMatType & ~MAT_COLL_PNT_OUTPUT);

          return GetIntegrationRule(eMatTypeColloc).iGetNumEvalPoints();
     }
protected:
     virtual const IntegrationRule& GetIntegrationRule(FemMatrixType eMatType) const=0;
     virtual Matrix InterpGaussToNodal(FemMatrixType eMatType, const Matrix& taun) const=0;
     virtual ComplexMatrix InterpGaussToNodal(FemMatrixType eMatType, const ComplexMatrix& taun) const=0;
     virtual void ScalarInterpMatrix(const ColumnVector& rv, Matrix& Hs, octave_idx_type irow) const=0;
     virtual void ScalarInterpMatrixDer(const ColumnVector& rv, Matrix& Hd) const {
          throw std::runtime_error("function ScalarInterpMatrixDer not implemented for this element type");
     }

     virtual void ScalarGradientMatrix(const ColumnVector& rv, const Matrix& J, double detJ, Matrix& invJ, Matrix& Bt) const {
          const octave_idx_type iNumNodes = X.columns();
          const octave_idx_type iNumDir = 3;
          FEM_ASSERT(rv.numel() == iNumDir);
          FEM_ASSERT(J.rows() == iNumDir);
          FEM_ASSERT(J.columns() == iNumDir);
          FEM_ASSERT(invJ.rows() == iNumDir);
          FEM_ASSERT(invJ.columns() == iNumDir);
          FEM_ASSERT(Bt.rows() == iNumDir);
          FEM_ASSERT(Bt.columns() == iNumNodes);

          Inverse3x3(J, detJ, invJ);

          Matrix Hd(iNumNodes, iNumDir);

          ScalarInterpMatrixDer(rv, Hd);

          for (octave_idx_type j = 0; j < iNumNodes; ++j) {
               for (octave_idx_type i = 0; i < iNumDir; ++i) {
                    double Bt_ij = 0.;

                    for (octave_idx_type k = 0; k < iNumDir; ++k) {
                         Bt_ij += invJ.xelem(i, k) * Hd.xelem(j, k);
                    }

                    Bt.xelem(i, j) = Bt_ij;
               }
          }
     }

     virtual void DispInterpMatrix(const ColumnVector& rv, Matrix& H) const {
          FEM_ASSERT(rv.numel() == 3);
          FEM_ASSERT(H.rows() == 3);
          FEM_ASSERT(H.columns() == nodes.numel() * 3);

          const octave_idx_type N = nodes.numel();

          Matrix Hs(1, N);

          ScalarInterpMatrix(rv, Hs, 0);

          std::fill(H.fortran_vec(), H.fortran_vec() + H.numel(), 0.);

          for (octave_idx_type i = 0; i < N; ++i) {
               for (octave_idx_type j = 0; j < 3; ++j) {
                    H.xelem(j, 3 * i + j) = Hs.xelem(i);
               }
          }
     }

     void ScalarInterpMatrixGrad(const ColumnVector& rv, const Matrix& invJ, const Matrix& Hd, Matrix& gradH) const {
          FEM_ASSERT(rv.rows() == 3);
          FEM_ASSERT(invJ.rows() == Hd.columns());
          FEM_ASSERT(gradH.rows() == Hd.rows());
          FEM_ASSERT(gradH.columns() == invJ.columns());
          FEM_ASSERT(Hd.rows() == X.columns());

          const octave_idx_type N = Hd.rows();
          const octave_idx_type M = invJ.columns();
          const octave_idx_type P = invJ.rows();

          for (octave_idx_type i = 0; i < N; ++i) {
               for (octave_idx_type j = 0; j < M; ++j) {
                    double gradHij = 0.;

                    for (octave_idx_type k = 0; k < P; ++k) {
                         gradHij += Hd.xelem(i, k) * invJ.xelem(j, k);
                    }

                    gradH.xelem(i, j) = gradHij;
               }
          }
     }

     void StrainMatrixNL(const ColumnVector& rv, const Matrix& gradH, Matrix& BNL) const {
          FEM_ASSERT(rv.rows() == 3);
          FEM_ASSERT(gradH.rows() == X.columns());
          FEM_ASSERT(gradH.columns() == 3);
          FEM_ASSERT(BNL.columns() == 3 * X.columns());
          FEM_ASSERT(BNL.rows() == 9);

          const octave_idx_type N = gradH.rows();

          BNL.fill(0.);

          for (octave_idx_type i = 0; i < 3; ++i) {
               for (octave_idx_type j = 0; j < 3; ++j) {
                    for (octave_idx_type k = 0; k < N; ++k) {
                         BNL.xelem(i * 3 + j, i + 3 * k) = gradH.xelem(k, j);
                    }
               }
          }
     }

     void StiffnessMatrixNL(Matrix& KNLe, MeshInfo& oMeshInfo, FemMatrixType eMatType) const {
          const IntegrationRule& oIntegRule = GetIntegrationRule(eMatType);
          const octave_idx_type iNumGauss = oIntegRule.iGetNumEvalPoints();
          const octave_idx_type iNumDof = iGetNumDof(eMatType);
          const octave_idx_type iNumDir = oIntegRule.iGetNumDirections();
          const octave_idx_type iNumNodes = nodes.numel();
          ColumnVector rv(iNumDir);
          constexpr octave_idx_type iNumStressComp = 6;
          constexpr octave_idx_type iNumStressCompNL = 9;

          FEM_ASSERT(id >= 1);
          FEM_ASSERT(iNumDof == 3 * X.columns());
          FEM_ASSERT(KNLe.rows() == iNumDof);
          FEM_ASSERT(KNLe.columns() == KNLe.rows());
          FEM_ASSERT(iNumDir == 3); // FIXME: Jacobian must be a 3x3 matrix

          Matrix J(iNumDir, iNumDir), invJ(iNumDir, iNumDir), H(1, iNumNodes);
          Matrix BNL(iNumStressCompNL, iNumDof), Hd(iNumNodes, iNumDir), gradH(iNumNodes, iNumDir); // FIXME: will not work for Tet10 elements!!
          Matrix SgBNL(iNumStressCompNL, iNumDof);
          ColumnVector taug(iNumStressComp);
          Matrix Sg(iNumStressCompNL, iNumStressCompNL);

          for (octave_idx_type i = 0; i < iNumGauss; ++i) {
               const double alpha = oIntegRule.dGetWeight(i);

               for (octave_idx_type j = 0; j < iNumDir; ++j) {
                    rv.xelem(j) = oIntegRule.dGetPosition(i, j);
               }

               const double detJ = Jacobian(rv, J);

               Inverse3x3(J, detJ, invJ);

               ScalarInterpMatrix(rv, H, 0);

               for (octave_idx_type j = 0; j < iNumStressComp; ++j) {
                    double taugj = 0.;

                    for (octave_idx_type k = 0; k < iNumNodes; ++k) {
                         taugj += H.xelem(k) * tauRef.xelem(j + iNumStressComp * k);
                    }

                    taug.xelem(j) = taugj;
               }

               for (octave_idx_type j = 0; j < iNumStressCompNL * iNumStressCompNL; ++j) {
                    Sg.xelem(j) = 0.;
               }

               for (octave_idx_type j = 0; j < 3; ++j) {
                    for (octave_idx_type k = 0; k < 3; ++k) {
                         Sg.xelem(j * 3 + k, j * 3 + k) = taug.xelem(k);
                    }

                    Sg.xelem(j * 3,     j * 3 + 1) = Sg.xelem(j * 3 + 1,     j * 3) = taug.xelem(3);
                    Sg.xelem(j * 3,     j * 3 + 2) = Sg.xelem(j * 3 + 2,     j * 3) = taug.xelem(5);
                    Sg.xelem(j * 3 + 1, j * 3 + 2) = Sg.xelem(j * 3 + 2, j * 3 + 1) = taug.xelem(4);
               }

               ScalarInterpMatrixDer(rv, Hd);
               ScalarInterpMatrixGrad(rv, invJ, Hd, gradH);
               StrainMatrixNL(rv, gradH, BNL);

               for (octave_idx_type l = 0; l < iNumStressCompNL; ++l) {
                    for (octave_idx_type m = 0; m < iNumDof; ++m) {
                         double SgBNLlm = 0.;

                         for (octave_idx_type n = 0; n < iNumStressCompNL; ++n) {
                              SgBNLlm += detJ * alpha * Sg.xelem(l + iNumStressCompNL * n) * BNL.xelem(n + iNumStressCompNL * m);
                         }

                         SgBNL.xelem(l + iNumStressCompNL * m) = SgBNLlm;
                    }
               }

               for (octave_idx_type l = 0; l < iNumDof; ++l) {
                    for (octave_idx_type m = l; m < iNumDof; ++m) {
                         double KNLelm = 0.;

                         for (octave_idx_type n = 0; n < iNumStressCompNL; ++n) {
                              KNLelm += BNL.xelem(n + iNumStressCompNL * l) * SgBNL.xelem(n + iNumStressCompNL * m);
                         }

                         KNLe.xelem(l + iNumDof * m) += KNLelm;
                    }
               }
          }

          for (octave_idx_type i = 1; i < iNumDof; ++i) {
               for (octave_idx_type j = 0; j < i; ++j) {
                    KNLe.xelem(i + iNumDof * j) = KNLe.xelem(j + iNumDof * i);
               }
          }
     }

     void AddMeshInfo(MeshInfo& info, const IntegrationRule& oIntegRule, double detJ) const {
          info.Add(MeshInfo::JACOBIAN_DET, detJ);
     }

     void StiffnessMatrix(Matrix& Ke, MeshInfo& info, FemMatrixType eMatType) const {
          const IntegrationRule& oIntegRule = GetIntegrationRule(eMatType);
          const octave_idx_type iNumDof = iGetNumDof(eMatType);
          const octave_idx_type iNumDir = oIntegRule.iGetNumDirections();
          ColumnVector rv(iNumDir);
          const Matrix& C = material->LinearElasticity();
          const octave_idx_type iNumStrains = C.rows();

          FEM_ASSERT(C.rows() == C.columns());
          FEM_ASSERT(Ke.rows() == iNumDof);
          FEM_ASSERT(Ke.columns() == iNumDof);

          Matrix J(iNumDir, iNumDir), invJ(iNumDir, iNumDir), B(iNumStrains, iNumDof), CB(iNumStrains, iNumDof);

          for (octave_idx_type i = 0; i < oIntegRule.iGetNumEvalPoints(); ++i) {
               const double alpha = oIntegRule.dGetWeight(i);

               for (octave_idx_type j = 0; j < iNumDir; ++j) {
                    rv.xelem(j) = oIntegRule.dGetPosition(i, j);
               }

               const double detJ = Jacobian(rv, J);

               AddMeshInfo(info, oIntegRule, detJ);

               StrainMatrix(rv, J, detJ, invJ, B);

               for (octave_idx_type l = 0; l < iNumStrains; ++l) {
                    for (octave_idx_type m = 0; m < iNumDof; ++m) {
                         double CBlm = 0.;

                         for (octave_idx_type n = 0; n < iNumStrains; ++n) {
                              CBlm += detJ * alpha * C.xelem(l + iNumStrains * n) * B.xelem(n + iNumStrains * m);
                         }

                         FEM_ASSERT(std::isfinite(CBlm));

                         CB.xelem(l + iNumStrains * m) = CBlm;
                    }
               }

#ifdef DEBUG
               for (octave_idx_type k = 0; k < CB.rows(); ++k) {
                    for (octave_idx_type j = 0; j < CB.columns(); ++j) {
                         FEM_ASSERT(std::isfinite(CB(k, j)));
                    }
               }
#endif

               for (octave_idx_type l = 0; l < iNumDof; ++l) {
                    for (octave_idx_type m = l; m < iNumDof; ++m) {
                         double Kelm = 0.;

                         for (octave_idx_type n = 0; n < iNumStrains; ++n) {
                              Kelm += B.xelem(n + iNumStrains * l) * CB.xelem(n + iNumStrains * m);
                         }

                         FEM_ASSERT(std::isfinite(Kelm));

                         Ke.xelem(l + iNumDof * m) += Kelm;
                    }
               }
          }

          for (octave_idx_type i = 1; i < iNumDof; ++i) {
               for (octave_idx_type j = 0; j < i; ++j) {
                    Ke.xelem(i + iNumDof * j) = Ke.xelem(j + iNumDof * i);
               }
          }

#ifdef DEBUG
          for (octave_idx_type i = 0; i < Ke.rows(); ++i) {
               for (octave_idx_type j = 0; j < Ke.columns(); ++j) {
                    FEM_ASSERT(std::isfinite(Ke(i, j)));
               }
          }
#endif
     }

     void MassMatrix(Matrix& Me, MeshInfo& info, FemMatrixType eMatType) const {
          const IntegrationRule& oIntegRule = GetIntegrationRule(eMatType);
          const octave_idx_type iNumDof = iGetNumDof(eMatType);
          const octave_idx_type iNumDir = oIntegRule.iGetNumDirections();
          ColumnVector rv(iNumDir);
          const double rho = material->Density();
          const octave_idx_type iNumDisp = X.rows();

          FEM_ASSERT(Me.rows() == iNumDof);
          FEM_ASSERT(Me.columns() == iNumDof);

          Matrix J(iNumDir, iNumDir), H(iNumDisp, iNumDof);

          for (octave_idx_type i = 0; i < oIntegRule.iGetNumEvalPoints(); ++i) {
               const double alpha = oIntegRule.dGetWeight(i);

               for (octave_idx_type j = 0; j < iNumDir; ++j) {
                    rv.xelem(j) = oIntegRule.dGetPosition(i, j);
               }

               const double detJ = Jacobian(rv, J);

               AddMeshInfo(info, oIntegRule, detJ);

               DispInterpMatrix(rv, H);

               for (octave_idx_type l = 0; l < iNumDof; ++l) {
                    for (octave_idx_type m = l; m < iNumDof; ++m) {
                         double Melm = 0.;

                         for (octave_idx_type n = 0; n < iNumDisp; ++n) {
                              Melm += H.xelem(n + iNumDisp * l) * H.xelem(n + iNumDisp * m);
                         }

                         FEM_ASSERT(std::isfinite(Melm));

                         Me.xelem(l + iNumDof * m) += Melm * alpha * rho * detJ;
                    }
               }
          }

          for (octave_idx_type i = 1; i < iNumDof; ++i) {
               for (octave_idx_type j = 0; j < i; ++j) {
                    Me.xelem(i + iNumDof * j) = Me.xelem(j + iNumDof * i);
               }
          }
#ifdef DEBUG
          for (octave_idx_type i = 0; i < Me.rows(); ++i) {
               for (octave_idx_type j = 0; j < Me.columns(); ++j) {
                    FEM_ASSERT(std::isfinite(Me(i, j)));
               }
          }
#endif
     }

     void DampingMatrix(Matrix& De, MeshInfo& info, FemMatrixType eMatType) const {
          const octave_idx_type iNumDof = iGetNumDof(eMatType);

          FEM_ASSERT(De.rows() == iNumDof);
          FEM_ASSERT(De.columns() == iNumDof);

          const double alpha = material->AlphaDamping();

          if (alpha) {
               MassMatrix(De, info, static_cast<FemMatrixType>(MAT_MASS | (eMatType & MAT_SYM_MASK)));

               for (octave_idx_type j = 0; j < iNumDof; ++j) {
                    for (octave_idx_type i = 0; i < iNumDof; ++i) {
                         De.xelem(i + iNumDof * j) *= alpha;
                    }
               }
          }

          const double beta = material->BetaDamping();

          if (beta) {
               Matrix Ke(iNumDof, iNumDof, 0.);

               StiffnessMatrix(Ke, info, static_cast<FemMatrixType>(MAT_STIFFNESS | (eMatType & MAT_SYM_MASK)));

               for (octave_idx_type j = 0; j < iNumDof; ++j) {
                    for (octave_idx_type i = 0; i < iNumDof; ++i) {
                         De.xelem(i + iNumDof * j) += beta * Ke.xelem(i + iNumDof * j);
                    }
               }
          }
     }

     void StiffnessMatrixTanDelta(Matrix& Ke_im, MeshInfo& info, FemMatrixType eMatType) const {
          const octave_idx_type iNumDof = iGetNumDof(eMatType);

          FEM_ASSERT(Ke_im.rows() == iNumDof);
          FEM_ASSERT(Ke_im.columns() == iNumDof);

          const double tan_delta = material->TanDeltaDamping();

          if (tan_delta) {
               StiffnessMatrix(Ke_im, info, static_cast<FemMatrixType>(MAT_STIFFNESS | (eMatType & MAT_SYM_MASK)));

               for (octave_idx_type j = 0; j < iNumDof; ++j) {
                    for (octave_idx_type i = 0; i < iNumDof; ++i) {
                         Ke_im.xelem(i + iNumDof * j) *= tan_delta;
                    }
               }
          }
     }

     void GravityLoadVector(Matrix& R, const Matrix& g, MeshInfo& info, FemMatrixType eMatType) const {
          const IntegrationRule& oIntegRule = GetIntegrationRule(eMatType);
          const octave_idx_type iNumDof = iGetNumDof(eMatType);
          const octave_idx_type iNumDir = oIntegRule.iGetNumDirections();
          ColumnVector rv(iNumDir);
          const double rho = material->Density();
          const octave_idx_type iNumDisp = X.rows();

          FEM_ASSERT(R.rows() == iNumDof);
          FEM_ASSERT(R.columns() == g.columns());
          FEM_ASSERT(g.rows() == iNumDisp);

          Matrix J(iNumDir, iNumDir), H(iNumDisp, iNumDof);

          for (octave_idx_type i = 0; i < oIntegRule.iGetNumEvalPoints(); ++i) {
               const double alpha = oIntegRule.dGetWeight(i);

               for (octave_idx_type j = 0; j < iNumDir; ++j) {
                    rv.xelem(j) = oIntegRule.dGetPosition(i, j);
               }

               const double detJ = Jacobian(rv, J);

               AddMeshInfo(info, oIntegRule, detJ);

               DispInterpMatrix(rv, H);

               const double dm = alpha * rho * detJ;

               for (octave_idx_type l = 0; l < g.columns(); ++l) {
                    for (octave_idx_type j = 0; j < iNumDof; ++j) {
                         double Rijl = 0.;

                         for (octave_idx_type k = 0; k < iNumDisp; ++k) {
                              Rijl += H.xelem(k + iNumDisp * j) * dm * g.xelem(k + iNumDisp * l);
                         }

                         R.xelem(j + iNumDof * l) += Rijl;
                    }
               }
          }
     }

     void CentrifugalLoadVector(Matrix& R, const NDArray& WxWx, MeshInfo& info, FemMatrixType eMatType) const {
          const IntegrationRule& oIntegRule = GetIntegrationRule(eMatType);
          const octave_idx_type iNumDof = iGetNumDof(eMatType);
          const octave_idx_type iNumDir = oIntegRule.iGetNumDirections();
          ColumnVector rv(iNumDir);
          const double rho = material->Density();
          const octave_idx_type iNumDisp = X.rows();
          const octave_idx_type iNumNodes = X.columns();
          const octave_idx_type iNumLoads = WxWx.ndims() >= 3 ? WxWx.dim3() : 1;
          FEM_ASSERT(R.rows() == iNumDof);
          FEM_ASSERT(R.columns() == iNumLoads);
          FEM_ASSERT(WxWx.rows() == iNumDisp);

          Matrix J(iNumDir, iNumDir), H(iNumDisp, iNumDof);
          ColumnVector Xi(iNumDisp), ai(iNumDisp);

          for (octave_idx_type i = 0; i < oIntegRule.iGetNumEvalPoints(); ++i) {
               const double alpha = oIntegRule.dGetWeight(i);

               for (octave_idx_type j = 0; j < iNumDir; ++j) {
                    rv.xelem(j) = oIntegRule.dGetPosition(i, j);
               }

               const double detJ = Jacobian(rv, J);

               AddMeshInfo(info, oIntegRule, detJ);

               DispInterpMatrix(rv, H);

               for (octave_idx_type k = 0; k < iNumDisp; ++k) {
                    double Xik = 0.;

                    for (octave_idx_type j = 0; j < iNumNodes; ++j) {
                         Xik += H.xelem(k, j * 3 + k) * X.xelem(k, j);
                    }

                    Xi.xelem(k) = Xik;
               }

               const double dm = alpha * rho * detJ;

               for (octave_idx_type l = 0; l < iNumLoads; ++l) {
                    for (octave_idx_type j = 0; j < iNumDisp; ++j) {
                         double aij = 0.;

                         for (octave_idx_type k = 0; k < iNumDisp; ++k) {
                              aij += WxWx.xelem(j, k, l) * Xi.xelem(k);
                         }

                         ai.xelem(j) = aij;
                    }

                    for (octave_idx_type j = 0; j < iNumDof; ++j) {
                         double Rijl = 0.;

                         for (octave_idx_type k = 0; k < iNumDisp; ++k) {
                              Rijl -= H.xelem(k + iNumDisp * j) * dm * ai.xelem(k);
                         }

                         R.xelem(j + iNumDof * l) += Rijl;
                    }
               }
          }
     }

     void CentrifugalStiffnessMatrix(Matrix& Ke, MeshInfo& info, FemMatrixType eMatType) const {
          const double rho = material->Density();

          RotatingRefFrameMassMatrix(oElemData.oCentrifugal.WxWx, rho, Ke, info, eMatType);
     }

     void CoriolisDampingMatrix(Matrix& De, MeshInfo& info, FemMatrixType eMatType) const {
          const double rho = material->Density();

          RotatingRefFrameMassMatrix(oElemData.oCentrifugal.Wx, 2. * rho, De, info, eMatType);
     }

     void AngularAccelerationStiffnessMatrix(Matrix& Ke, MeshInfo& info, FemMatrixType eMatType) const {
          const double rho = material->Density();

          RotatingRefFrameMassMatrix(oElemData.oCentrifugal.WPx, rho, Ke, info, eMatType);
     }

     void RotatingRefFrameMassMatrix(const NDArray& AW, const double beta, Matrix& Ke, MeshInfo& info, FemMatrixType eMatType) const {
          const IntegrationRule& oIntegRule = GetIntegrationRule(eMatType);
          const octave_idx_type iNumDof = iGetNumDof(eMatType);
          const octave_idx_type iNumDir = oIntegRule.iGetNumDirections();
          ColumnVector rv(iNumDir);
          const octave_idx_type iNumDisp = X.rows();

          FEM_ASSERT(iNumDisp == AW.rows());
          FEM_ASSERT(iNumDisp == AW.columns());
          FEM_ASSERT(Ke.rows() == iNumDof);
          FEM_ASSERT(Ke.columns() == iNumDof);

          Matrix J(iNumDir, iNumDir), H(iNumDisp, iNumDof), AWH(iNumDisp, iNumDof);

          for (octave_idx_type i = 0; i < oIntegRule.iGetNumEvalPoints(); ++i) {
               const double alpha = oIntegRule.dGetWeight(i);

               for (octave_idx_type j = 0; j < iNumDir; ++j) {
                    rv.xelem(j) = oIntegRule.dGetPosition(i, j);
               }

               const double detJ = Jacobian(rv, J);

               AddMeshInfo(info, oIntegRule, detJ);

               DispInterpMatrix(rv, H);

               for (octave_idx_type l = 0; l < iNumDof; ++l) {
                    for (octave_idx_type j = 0; j < iNumDisp; ++j) {
                         double AWHjl = 0.;

                         for (octave_idx_type k = 0; k < iNumDisp; ++k) {
                              AWHjl += AW.xelem(j, k) * H.xelem(k, l);
                         }

                         AWH.xelem(j, l) = AWHjl;
                    }
               }

               for (octave_idx_type l = 0; l < iNumDof; ++l) {
                    for (octave_idx_type m = 0; m < iNumDof; ++m) {
                         double Kelm = 0.;

                         for (octave_idx_type n = 0; n < iNumDisp; ++n) {
                              Kelm += H.xelem(n + iNumDisp * l) * AWH.xelem(n + iNumDisp * m);
                         }

                         Ke.xelem(l + iNumDof * m) += Kelm * alpha * beta * detJ;
                    }
               }
          }

#ifdef DEBUG
          for (octave_idx_type i = 0; i < Ke.rows(); ++i) {
               for (octave_idx_type j = 0; j < Ke.columns(); ++j) {
                    FEM_ASSERT(std::isfinite(Ke(i, j)));
               }
          }
#endif
     }

     void StructuralLoadVector(Matrix& R, MeshInfo& info, FemMatrixType eMatType) const {
          if (oElemData.oGravity.g.columns()) {
               GravityLoadVector(R, oElemData.oGravity.g, info, eMatType);
          }

          if (oElemData.oCentrifugal.WxWx.numel()) {
               CentrifugalLoadVector(R, oElemData.oCentrifugal.WxWx, info, eMatType);
          }

          if (oElemData.oCentrifugal.WPx.numel()) {
               CentrifugalLoadVector(R, oElemData.oCentrifugal.WPx, info, eMatType);
          }

          if (epsilonRef.columns() || dTheta.columns()) {
               ThermalLoadVector(R, info, eMatType);
          }
     }

     void InertiaMoment1(NDArray& S, FemMatrixType eMatType) const {
          const IntegrationRule& oIntegRule = GetIntegrationRule(eMatType);
          const octave_idx_type iNumDof = iGetNumDof(eMatType);
          const octave_idx_type iNumDir = oIntegRule.iGetNumDirections();
          const double rho = material->Density();
          const octave_idx_type iNumDisp = X.rows();
          const octave_idx_type iNumNodes = nodes.numel();
          ColumnVector rv(iNumDir);

          constexpr octave_idx_type Srows = 3;

          FEM_ASSERT(X.rows() == 3);
          FEM_ASSERT(S.ndims() == 2);
          FEM_ASSERT(S.rows() == Srows);
          FEM_ASSERT(S.columns() == 1);
          FEM_ASSERT(iNumDisp == S.rows());

          Matrix J(iNumDir, iNumDir), H(iNumDisp, iNumDof);

          for (octave_idx_type i = 0; i < oIntegRule.iGetNumEvalPoints(); ++i) {
               const double alpha = oIntegRule.dGetWeight(i);

               for (octave_idx_type j = 0; j < iNumDir; ++j) {
                    rv.xelem(j) = oIntegRule.dGetPosition(i, j);
               }

               const double detJ = Jacobian(rv, J);

               DispInterpMatrix(rv, H);

               for (octave_idx_type l = 0; l < Srows; ++l) {
                    double fil = 0.;

                    for (octave_idx_type n = 0; n < iNumNodes; ++n) {
                         for (octave_idx_type m = 0; m < iNumDisp; ++m) {
                              fil += H.xelem(l + iNumDisp * (iNumDisp * n + m)) * X.xelem(m + 3 * n);
                         }
                    }

                    atomic_fetch_add_float(S.xelem(l), fil * alpha * rho * detJ);
               }
          }
     }

     void InertiaMatrix(NDArray& Inv7, FemMatrixType eMatType) const {
          const IntegrationRule& oIntegRule = GetIntegrationRule(eMatType);
          const octave_idx_type iNumDof = iGetNumDof(eMatType);
          const octave_idx_type iNumDir = oIntegRule.iGetNumDirections();
          const double rho = material->Density();
          const octave_idx_type iNumDisp = X.rows();
          const octave_idx_type iNumNodes = nodes.numel();
          ColumnVector rv(iNumDir);

          FEM_ASSERT(X.rows() == iNumDisp);
          FEM_ASSERT(iNumDisp == 3);
          FEM_ASSERT(Inv7.ndims() == 2);
          FEM_ASSERT(Inv7.rows() == 3);
          FEM_ASSERT(Inv7.columns() == 3);
          FEM_ASSERT(iNumDisp == Inv7.rows());

          Matrix J(iNumDir, iNumDir), H(iNumDisp, iNumDof);
          ColumnVector fi(iNumDisp);

          for (octave_idx_type i = 0; i < oIntegRule.iGetNumEvalPoints(); ++i) {
               const double alpha = oIntegRule.dGetWeight(i);

               for (octave_idx_type j = 0; j < iNumDir; ++j) {
                    rv.xelem(j) = oIntegRule.dGetPosition(i, j);
               }

               const double detJ = Jacobian(rv, J);

               DispInterpMatrix(rv, H);

               for (octave_idx_type l = 0; l < iNumDisp; ++l) {
                    fi.xelem(l) = 0.;

                    for (octave_idx_type n = 0; n < iNumNodes; ++n) {
                         for (octave_idx_type m = 0; m < iNumDisp; ++m) {
                              fi.xelem(l) += H.xelem(l + iNumDisp * (iNumDisp * n + m)) * X.xelem(m + iNumDisp * n);
                         }
                    }
               }

               const double dmi = alpha * rho * detJ;

               atomic_fetch_add_float(Inv7.xelem(0 + 3 * 0), (fi.xelem(1) * fi.xelem(1) + fi.xelem(2) * fi.xelem(2)) * dmi);
               atomic_fetch_add_float(Inv7.xelem(0 + 3 * 1), -(fi.xelem(0) * fi.xelem(1)) * dmi);
               atomic_fetch_add_float(Inv7.xelem(0 + 3 * 2), -(fi.xelem(0) * fi.xelem(2)) * dmi);
               atomic_fetch_add_float(Inv7.xelem(1 + 3 * 1), (fi.xelem(0) * fi.xelem(0) + fi.xelem(2) * fi.xelem(2)) * dmi);
               atomic_fetch_add_float(Inv7.xelem(1 + 3 * 2), -(fi.xelem(1) * fi.xelem(2)) * dmi);
               atomic_fetch_add_float(Inv7.xelem(2 + 3 * 2), (fi.xelem(0) * fi.xelem(0) + fi.xelem(1) * fi.xelem(1)) * dmi);
          }
     }

     void InertiaInv3(NDArray& Inv3, FemMatrixType eMatType, const NDArray& U) const {
          const IntegrationRule& oIntegRule = GetIntegrationRule(eMatType);
          const octave_idx_type iNumDof = iGetNumDof(eMatType);
          const octave_idx_type iNumDir = oIntegRule.iGetNumDirections();
          const double rho = material->Density();
          const octave_idx_type iNumDisp = X.rows();
          const octave_idx_type iNumNodes = nodes.numel();
          const octave_idx_type iNumNodesMesh = U.dim1();
          const octave_idx_type iNumModes = U.ndims() >= 3 ? U.dim3() : 1;

          ColumnVector rv(iNumDir);

          FEM_ASSERT(U.ndims() == 3);
          FEM_ASSERT(U.dim2() >= iNumDisp);
          FEM_ASSERT(U.dim2() == 6);
          FEM_ASSERT(Inv3.ndims() == 2);
          FEM_ASSERT(Inv3.rows() == 3);
          FEM_ASSERT(Inv3.columns() == iNumModes);
          FEM_ASSERT(iNumDisp == Inv3.rows());

          Matrix J(iNumDir, iNumDir), H(iNumDisp, iNumDof);

          for (octave_idx_type i = 0; i < oIntegRule.iGetNumEvalPoints(); ++i) {
               const double alpha = oIntegRule.dGetWeight(i);

               for (octave_idx_type j = 0; j < iNumDir; ++j) {
                    rv.xelem(j) = oIntegRule.dGetPosition(i, j);
               }

               const double detJ = Jacobian(rv, J);

               DispInterpMatrix(rv, H);

               const double dmi = alpha * rho * detJ;

               for (octave_idx_type j = 0; j < iNumModes; ++j) {
                    for (octave_idx_type l = 0; l < iNumDisp; ++l) {
                         double Uil = 0.;

                         for (octave_idx_type m = 0; m < iNumDisp; ++m) {
                              for (octave_idx_type n = 0; n < iNumNodes; ++n) {
                                   Uil += H.xelem(l + iNumDisp * (iNumDisp * n + m)) * U.xelem(nodes.xelem(n).value() - 1 + iNumNodesMesh * (m + 6 * j));
                              }
                         }

                         atomic_fetch_add_float(Inv3.xelem(l + 3 * j), dmi * Uil);
                    }
               }
          }
     }

     void InertiaInv4(NDArray& Inv4, FemMatrixType eMatType, const NDArray& U) const {
          const IntegrationRule& oIntegRule = GetIntegrationRule(eMatType);
          const octave_idx_type iNumDof = iGetNumDof(eMatType);
          const octave_idx_type iNumDir = oIntegRule.iGetNumDirections();
          const double rho = material->Density();
          const octave_idx_type iNumDisp = X.rows();
          const octave_idx_type iNumNodes = nodes.numel();
          const octave_idx_type iNumNodesMesh = U.dim1();
          const octave_idx_type iNumModes = U.ndims() >= 3 ? U.dim3() : 1;
          ColumnVector rv(iNumDir);

          FEM_ASSERT(U.ndims() == 3);
          FEM_ASSERT(U.dim2() >= iNumDisp);
          FEM_ASSERT(U.dim2() == 6);
          FEM_ASSERT(Inv4.rows() == 3);
          FEM_ASSERT(Inv4.columns() == iNumModes);
          FEM_ASSERT(iNumDisp == Inv4.rows());

          Matrix J(iNumDir, iNumDir), H(iNumDisp, iNumDof);
          ColumnVector Ui(iNumDisp), fi(iNumDisp);

          for (octave_idx_type i = 0; i < oIntegRule.iGetNumEvalPoints(); ++i) {
               const double alpha = oIntegRule.dGetWeight(i);

               for (octave_idx_type j = 0; j < iNumDir; ++j) {
                    rv.xelem(j) = oIntegRule.dGetPosition(i, j);
               }

               const double detJ = Jacobian(rv, J);

               DispInterpMatrix(rv, H);

               const double dmi = alpha * rho * detJ;

               for (octave_idx_type j = 0; j < iNumModes; ++j) {
                    for (octave_idx_type l = 0; l < iNumDisp; ++l) {
                         Ui.xelem(l) = 0.;
                         fi.xelem(l) = 0.;

                         for (octave_idx_type m = 0; m < iNumDisp; ++m) {
                              for (octave_idx_type n = 0; n < iNumNodes; ++n) {
                                   const double Hlmn = H.xelem(l + iNumDisp * (iNumDisp * n + m));
                                   Ui.xelem(l) += Hlmn * U.xelem(nodes.xelem(n).value() - 1 + iNumNodesMesh * (m + 6 * j));
                                   fi.xelem(l) += Hlmn * X.xelem(m + iNumDisp * n);
                              }
                         }
                    }

                    atomic_fetch_add_float(Inv4.xelem(0 + 3 * j), (fi.xelem(1) * Ui.xelem(2) - Ui.xelem(1) * fi.xelem(2)) * dmi);
                    atomic_fetch_add_float(Inv4.xelem(1 + 3 * j), (Ui.xelem(0) * fi.xelem(2) - fi.xelem(0) * Ui.xelem(2)) * dmi);
                    atomic_fetch_add_float(Inv4.xelem(2 + 3 * j), (fi.xelem(0) * Ui.xelem(1) - Ui.xelem(0) * fi.xelem(1)) * dmi);
               }
          }
     }

     void InertiaInv5(NDArray& Inv5, FemMatrixType eMatType, const NDArray& U) const {
          const IntegrationRule& oIntegRule = GetIntegrationRule(eMatType);
          const octave_idx_type iNumDof = iGetNumDof(eMatType);
          const octave_idx_type iNumDir = oIntegRule.iGetNumDirections();
          const double rho = material->Density();
          const octave_idx_type iNumDisp = X.rows();
          const octave_idx_type iNumNodes = nodes.numel();
          const octave_idx_type iNumGauss = oIntegRule.iGetNumEvalPoints();
          const octave_idx_type iNumModes = U.ndims() >= 3 ? U.dim3() : 1;
          const octave_idx_type iNumNodesMesh = U.dim1();

          ColumnVector rv(iNumDir);

          FEM_ASSERT(U.ndims() == 3);
          FEM_ASSERT(U.dim2() >= iNumDisp);
          FEM_ASSERT(U.dim2() == 6);
          FEM_ASSERT(Inv5.ndims() == 3);
          FEM_ASSERT(Inv5.dim1() == 3);
          FEM_ASSERT(Inv5.dim2() == iNumModes);
          FEM_ASSERT(Inv5.dim3() == iNumModes);
          FEM_ASSERT(iNumDisp == Inv5.dim1());

          Matrix J(iNumDir, iNumDir), H(iNumDisp, iNumDof);
          NDArray Ui(dim_vector(iNumDisp, iNumModes, iNumGauss), 0.);
          ColumnVector dmi(iNumGauss);

          for (octave_idx_type i = 0; i < iNumGauss; ++i) {
               const double alpha = oIntegRule.dGetWeight(i);

               for (octave_idx_type j = 0; j < iNumDir; ++j) {
                    rv.xelem(j) = oIntegRule.dGetPosition(i, j);
               }

               const double detJ = Jacobian(rv, J);

               DispInterpMatrix(rv, H);

               dmi.xelem(i) = alpha * rho * detJ;

               for (octave_idx_type j = 0; j < iNumModes; ++j) {
                    for (octave_idx_type l = 0; l < iNumDisp; ++l) {
                         for (octave_idx_type m = 0; m < iNumDisp; ++m) {
                              for (octave_idx_type n = 0; n < iNumNodes; ++n) {
                                   Ui.xelem(l + iNumDisp * (j + iNumModes * i)) += H.xelem(l + iNumDisp * (iNumDisp * n + m)) * U.xelem(nodes.xelem(n).value() - 1 + iNumNodesMesh * (m + 6 * j));
                              }
                         }
                    }
               }
          }

          for (octave_idx_type i = 0; i < iNumGauss; ++i) {
               for (octave_idx_type j = 0; j < iNumModes; ++j) {
                    for (octave_idx_type k = 0; k < iNumModes; ++k) {
                         atomic_fetch_add_float(Inv5.xelem(0 + 3 * (k + iNumModes * j)), dmi.xelem(i) * (Ui.xelem(1 + iNumDisp * (j + iNumModes * i)) * Ui.xelem(2 + iNumDisp * (k + iNumModes * i)) - Ui.xelem(2 + iNumDisp * (j + iNumModes * i)) * Ui.xelem(1 + iNumDisp * (k + iNumModes * i))));
                         atomic_fetch_add_float(Inv5.xelem(1 + 3 * (k + iNumModes * j)), dmi.xelem(i) * (Ui.xelem(2 + iNumDisp * (j + iNumModes * i)) * Ui.xelem(0 + iNumDisp * (k + iNumModes * i)) - Ui.xelem(0 + iNumDisp * (j + iNumModes * i)) * Ui.xelem(2 + iNumDisp * (k + iNumModes * i))));
                         atomic_fetch_add_float(Inv5.xelem(2 + 3 * (k + iNumModes * j)), dmi.xelem(i) * (Ui.xelem(0 + iNumDisp * (j + iNumModes * i)) * Ui.xelem(1 + iNumDisp * (k + iNumModes * i)) - Ui.xelem(1 + iNumDisp * (j + iNumModes * i)) * Ui.xelem(0 + iNumDisp * (k + iNumModes * i))));
                    }
               }
          }
     }

     void InertiaInv8(NDArray& Inv8, FemMatrixType eMatType, const NDArray& U) const {
          const IntegrationRule& oIntegRule = GetIntegrationRule(eMatType);
          const octave_idx_type iNumDof = iGetNumDof(eMatType);
          const octave_idx_type iNumDir = oIntegRule.iGetNumDirections();
          const double rho = material->Density();
          const octave_idx_type iNumDisp = X.rows();
          const octave_idx_type iNumNodes = nodes.numel();
          const octave_idx_type iNumNodesMesh = U.dim1();
          const octave_idx_type iNumModes = U.ndims() >= 3 ? U.dim3() : 1;
          ColumnVector rv(iNumDir);

          FEM_ASSERT(X.rows() == 3);
          FEM_ASSERT(U.ndims() == 3);
          FEM_ASSERT(U.dim2() >= iNumDisp);
          FEM_ASSERT(U.dim2() == 6);
          FEM_ASSERT(Inv8.ndims() == 3);
          FEM_ASSERT(Inv8.dim1() == 3);
          FEM_ASSERT(Inv8.dim2() == 3);
          FEM_ASSERT(Inv8.dim3() == iNumModes);
          FEM_ASSERT(iNumDisp == Inv8.dim1());

          Matrix J(iNumDir, iNumDir), H(iNumDisp, iNumDof);
          ColumnVector Ui(iNumDisp), fi(iNumDisp);

          for (octave_idx_type i = 0; i < oIntegRule.iGetNumEvalPoints(); ++i) {
               const double alpha = oIntegRule.dGetWeight(i);

               for (octave_idx_type j = 0; j < iNumDir; ++j) {
                    rv.xelem(j) = oIntegRule.dGetPosition(i, j);
               }

               const double detJ = Jacobian(rv, J);

               DispInterpMatrix(rv, H);

               const double dmi = alpha * rho * detJ;

               for (octave_idx_type j = 0; j < iNumModes; ++j) {
                    for (octave_idx_type l = 0; l < iNumDisp; ++l) {
                         Ui.xelem(l) = 0.;
                         fi.xelem(l) = 0.;

                         for (octave_idx_type m = 0; m < iNumDisp; ++m) {
                              for (octave_idx_type n = 0; n < iNumNodes; ++n) {
                                   const double Hlmn = H.xelem(l + iNumDisp * (iNumDisp * n + m));
                                   Ui.xelem(l) += Hlmn * U.xelem(nodes.xelem(n).value() - 1 + iNumNodesMesh * (m + 6 * j));
                                   fi.xelem(l) += Hlmn * X.xelem(m + iNumDisp * n);
                              }
                         }
                    }

                    atomic_fetch_add_float(Inv8.xelem(0 + 3 * (0 + 3 * j)), (fi.xelem(2) * Ui.xelem(2) + fi.xelem(1) * Ui.xelem(1)) * dmi);
                    atomic_fetch_add_float(Inv8.xelem(1 + 3 * (0 + 3 * j)), -fi.xelem(0) * Ui.xelem(1) * dmi);
                    atomic_fetch_add_float(Inv8.xelem(2 + 3 * (0 + 3 * j)), -fi.xelem(0) * Ui.xelem(2) * dmi);
                    atomic_fetch_add_float(Inv8.xelem(0 + 3 * (1 + 3 * j)), -Ui.xelem(0) * fi.xelem(1) * dmi);
                    atomic_fetch_add_float(Inv8.xelem(1 + 3 * (1 + 3 * j)), (fi.xelem(2) * Ui.xelem(2) + fi.xelem(0) * Ui.xelem(0)) * dmi);
                    atomic_fetch_add_float(Inv8.xelem(2 + 3 * (1 + 3 * j)), -fi.xelem(1) * Ui.xelem(2) * dmi);
                    atomic_fetch_add_float(Inv8.xelem(0 + 3 * (2 + 3 * j)), -Ui.xelem(0) * fi.xelem(2) * dmi);
                    atomic_fetch_add_float(Inv8.xelem(1 + 3 * (2 + 3 * j)), -Ui.xelem(1) * fi.xelem(2) * dmi);
                    atomic_fetch_add_float(Inv8.xelem(2 + 3 * (2 + 3 * j)), (fi.xelem(1) * Ui.xelem(1) + fi.xelem(0) * Ui.xelem(0)) * dmi);
               }
          }
     }

     void InertiaInv9(NDArray& Inv9, FemMatrixType eMatType, const NDArray& U) const {
          const IntegrationRule& oIntegRule = GetIntegrationRule(eMatType);
          const octave_idx_type iNumDof = iGetNumDof(eMatType);
          const octave_idx_type iNumDir = oIntegRule.iGetNumDirections();
          const double rho = material->Density();
          const octave_idx_type iNumDisp = X.rows();
          const octave_idx_type iNumNodes = nodes.numel();
          const octave_idx_type iNumNodesMesh = U.dim1();
          const octave_idx_type iNumGauss = oIntegRule.iGetNumEvalPoints();
          const octave_idx_type iNumModes = U.ndims() >= 3 ? U.dim3() : 1;
          ColumnVector rv(iNumDir);

          FEM_ASSERT(U.ndims() == 3);
          FEM_ASSERT(U.dim2() >= iNumDisp);
          FEM_ASSERT(U.dim2() == 6);
          FEM_ASSERT(Inv9.ndims() == 4);
          FEM_ASSERT(Inv9.dim1() == 3);
          FEM_ASSERT(Inv9.dim2() == 3);
          FEM_ASSERT(Inv9.dim3() == iNumModes);
          FEM_ASSERT(Inv9.dims()(3) == iNumModes);
          FEM_ASSERT(iNumDisp == Inv9.dim1());

          Matrix J(iNumDir, iNumDir), H(iNumDisp, iNumDof);
          NDArray Ui(dim_vector(iNumDisp, iNumModes, iNumGauss), 0.);
          ColumnVector dmi(iNumGauss);

          for (octave_idx_type i = 0; i < iNumGauss; ++i) {
               const double alpha = oIntegRule.dGetWeight(i);

               for (octave_idx_type j = 0; j < iNumDir; ++j) {
                    rv.xelem(j) = oIntegRule.dGetPosition(i, j);
               }

               const double detJ = Jacobian(rv, J);

               DispInterpMatrix(rv, H);

               dmi.xelem(i) = alpha * rho * detJ;

               for (octave_idx_type j = 0; j < iNumModes; ++j) {
                    for (octave_idx_type l = 0; l < iNumDisp; ++l) {
                         for (octave_idx_type m = 0; m < iNumDisp; ++m) {
                              for (octave_idx_type n = 0; n < iNumNodes; ++n) {
                                   Ui.xelem(l + iNumDisp * (j + iNumModes * i)) += H.xelem(l + iNumDisp * (iNumDisp * n + m)) * U.xelem(nodes.xelem(n).value() - 1 + iNumNodesMesh * (m + 6 * j));
                              }
                         }
                    }
               }
          }

          for (octave_idx_type i = 0; i < iNumGauss; ++i) {
               for (octave_idx_type j = 0; j < iNumModes; ++j) {
                    for (octave_idx_type k = 0; k < iNumModes; ++k) {
                         const octave_idx_type jk = j + k * Inv9.dim3();

                         atomic_fetch_add_float(Inv9.xelem(0 + 3 * (0 + 3 * jk)), (-Ui.xelem(2 + iNumDisp * (j + iNumModes * i)) * Ui.xelem(2 + iNumDisp * (k + iNumModes * i)) - Ui.xelem(1 + iNumDisp * (j + iNumModes * i)) * Ui.xelem(1 + iNumDisp * (k + iNumModes * i))) * dmi.xelem(i));
                         atomic_fetch_add_float(Inv9.xelem(1 + 3 * (0 + 3 * jk)), Ui.xelem(0 + iNumDisp * (j + iNumModes * i)) * Ui.xelem(1 + iNumDisp * (k + iNumModes * i)) * dmi.xelem(i));
                         atomic_fetch_add_float(Inv9.xelem(2 + 3 * (0 + 3 * jk)), Ui.xelem(0 + iNumDisp * (j + iNumModes * i)) * Ui.xelem(2 + iNumDisp * (k + iNumModes * i)) * dmi.xelem(i));
                         atomic_fetch_add_float(Inv9.xelem(0 + 3 * (1 + 3 * jk)), Ui.xelem(0 + iNumDisp * (k + iNumModes * i)) * Ui.xelem(1 + iNumDisp * (j + iNumModes * i)) * dmi.xelem(i));
                         atomic_fetch_add_float(Inv9.xelem(1 + 3 * (1 + 3 * jk)), (-Ui.xelem(2 + iNumDisp * (j + iNumModes * i)) * Ui.xelem(2 + iNumDisp * (k + iNumModes * i)) - Ui.xelem(0 + iNumDisp * (j + iNumModes * i)) * Ui.xelem(0 + iNumDisp * (k + iNumModes * i))) * dmi.xelem(i));
                         atomic_fetch_add_float(Inv9.xelem(2 + 3 * (1 + 3 * jk)), Ui.xelem(1 + iNumDisp * (j + iNumModes * i)) * Ui.xelem(2 + iNumDisp * (k + iNumModes * i)) * dmi.xelem(i));
                         atomic_fetch_add_float(Inv9.xelem(0 + 3 * (2 + 3 * jk)), Ui.xelem(0 + iNumDisp * (k + iNumModes * i)) * Ui.xelem(2 + iNumDisp * (j + iNumModes * i)) * dmi.xelem(i));
                         atomic_fetch_add_float(Inv9.xelem(1 + 3 * (2 + 3 * jk)), Ui.xelem(1 + iNumDisp * (k + iNumModes * i)) * Ui.xelem(2 + iNumDisp * (j + iNumModes * i)) * dmi.xelem(i));
                         atomic_fetch_add_float(Inv9.xelem(2 + 3 * (2 + 3 * jk)), (-Ui.xelem(1 + iNumDisp * (j + iNumModes * i)) * Ui.xelem(1 + iNumDisp * (k + iNumModes * i)) - Ui.xelem(0 + iNumDisp * (j + iNumModes * i)) * Ui.xelem(0 + iNumDisp * (k + iNumModes * i))) * dmi.xelem(i));
                    }
               }
          }
     }

     void StrainNodalElem(NDArray& epsilonn, FemMatrixType eMatType, const NDArray& U) const {
          const IntegrationRule& oIntegRule = GetIntegrationRule(eMatType);
          const octave_idx_type iNumGauss = oIntegRule.iGetNumEvalPoints();
          const octave_idx_type iNumDof = iGetNumDof(eMatType);
          const octave_idx_type iNumDir = oIntegRule.iGetNumDirections();
          const octave_idx_type iNumLoads = U.ndims() >= 3 ? U.dim3() : 1;
          const octave_idx_type iNumNodes = nodes.numel();
          const octave_idx_type iNumNodesMesh = U.dim1();
          const octave_idx_type iNumElem = epsilonn.dim1();

          ColumnVector rv(iNumDir);
          const octave_idx_type iNumStrains = material->LinearElasticity().rows();

          FEM_ASSERT(epsilonn.ndims() >= 3);
          FEM_ASSERT(id >= 1);
          FEM_ASSERT(epsilonn.dim1() >= id);
          FEM_ASSERT(epsilonn.dim2() == iNumNodes);
          FEM_ASSERT(epsilonn.dim3() == iNumStrains);
          FEM_ASSERT((iNumLoads == 1 && epsilonn.ndims() == 3) || (epsilonn.dims()(3) == iNumLoads));
          FEM_ASSERT(U.ndims() >= 2);
          FEM_ASSERT(U.dim2() >= 3);
          FEM_ASSERT(U.dim2() == 6);
          FEM_ASSERT(iNumDof == 3 * X.columns());

          Matrix J(iNumDir, iNumDir), invJ(iNumDir, iNumDir), B(iNumStrains, iNumDof);
          ColumnVector Ue(iNumDof);
          Matrix epsilong(iNumGauss, iNumStrains * iNumLoads);

          for (octave_idx_type i = 0; i < iNumGauss; ++i) {
               for (octave_idx_type j = 0; j < iNumDir; ++j) {
                    rv.xelem(j) = oIntegRule.dGetPosition(i, j);
               }

               const double detJ = Jacobian(rv, J);

               StrainMatrix(rv, J, detJ, invJ, B);

               for (octave_idx_type l = 0; l < iNumLoads; ++l) {
                    for (octave_idx_type j = 0; j < iNumNodes; ++j) {
                         for (octave_idx_type k = 0; k < 3; ++k) {
                              FEM_ASSERT(nodes(j).value() > 0);
                              FEM_ASSERT(nodes(j).value() <= U.dim1());
                              Ue.xelem(3 * j + k) = U.xelem(nodes.xelem(j).value() - 1 + iNumNodesMesh * (k + 6 * l));
                         }
                    }

                    for (octave_idx_type j = 0; j < iNumStrains; ++j) {
                         double epsilongj = 0.;

                         for (octave_idx_type k = 0; k < iNumDof; ++k) {
                              epsilongj += B.xelem(j + iNumStrains * k) * Ue.xelem(k);
                         }

                         epsilong.xelem(i + iNumGauss * (l * iNumStrains + j)) = epsilongj;
                    }
               }
          }

          const Matrix epsilonen = InterpGaussToNodal(eMatType, epsilong);

          FEM_ASSERT(epsilonen.rows() == iNumNodes);
          FEM_ASSERT(epsilonen.columns() == epsilong.columns());

          for (octave_idx_type k = 0; k < iNumLoads; ++k) {
               for (octave_idx_type j = 0; j < iNumStrains; ++j) {
                    for (octave_idx_type i = 0; i < iNumNodes; ++i) {
                         epsilonn.xelem(id - 1 + iNumElem * (i + iNumNodes * (j + k * iNumStrains))) = epsilonen.xelem(i + iNumNodes * (k * iNumStrains + j));
                    }
               }
          }
     }

     void StressNodalElem(NDArray& taun, FemMatrixType eMatType, const NDArray& U) const {
          const IntegrationRule& oIntegRule = GetIntegrationRule(eMatType);
          const octave_idx_type iNumGauss = oIntegRule.iGetNumEvalPoints();
          const octave_idx_type iNumDof = iGetNumDof(eMatType);
          const octave_idx_type iNumDir = oIntegRule.iGetNumDirections();
          const octave_idx_type iNumLoads = U.ndims() >= 3 ? U.dim3() : 1;
          const octave_idx_type iNumNodes = nodes.numel();
          const octave_idx_type iNumElem = taun.dim1();
          const octave_idx_type iNumNodesMesh = U.dim1();
          ColumnVector rv(iNumDir);
          const Matrix& C = material->LinearElasticity();
          const double gamma = material->ThermalExpansion();
          const octave_idx_type iNumStrains = C.rows();

          FEM_ASSERT(C.rows() == C.columns());
          FEM_ASSERT(taun.ndims() >= 3);
          FEM_ASSERT(id >= 1);
          FEM_ASSERT(taun.dim1() >= id);
          FEM_ASSERT(taun.dim2() == iNumNodes);
          FEM_ASSERT(taun.dim3() == iNumStrains);
          FEM_ASSERT((iNumLoads == 1 && taun.ndims() == 3) || (taun.dims()(3) == iNumLoads));
          FEM_ASSERT(U.ndims() >= 2);
          FEM_ASSERT(U.dim2() >= 3);
          FEM_ASSERT(U.dim2() == 6);
          FEM_ASSERT(iNumDof == 3 * X.columns());

          Matrix J(iNumDir, iNumDir), invJ(iNumDir, iNumDir), B(iNumStrains, iNumDof), CB(iNumStrains, iNumDof), Ht;
          ColumnVector Ue(iNumDof), epsilonik(iNumPreLoads ? iNumStrains : 0);
          Matrix taug(iNumGauss, iNumStrains * iNumLoads);

          if (iNumPreLoads) {
               Ht.resize(1, iNumNodes);
          }

          for (octave_idx_type i = 0; i < iNumGauss; ++i) {
               for (octave_idx_type j = 0; j < iNumDir; ++j) {
                    rv.xelem(j) = oIntegRule.dGetPosition(i, j);
               }

               const double detJ = Jacobian(rv, J);

               StrainMatrix(rv, J, detJ, invJ, B);

               if (iNumPreLoads) {
                    ScalarInterpMatrix(rv, Ht, 0);
               }

               for (octave_idx_type j = 0; j < iNumStrains; ++j) {
                    for (octave_idx_type k = 0; k < iNumDof; ++k) {
                         double CBjk = 0;

                         for (octave_idx_type l = 0; l < iNumStrains; ++l) {
                              CBjk += C.xelem(j + iNumStrains * l) * B.xelem(l + iNumStrains * k);
                         }

                         CB.xelem(j + iNumStrains * k) = CBjk;
                    }
               }

               for (octave_idx_type l = 0; l < iNumLoads; ++l) {
                    for (octave_idx_type j = 0; j < iNumNodes; ++j) {
                         for (octave_idx_type k = 0; k < 3; ++k) {
                              FEM_ASSERT(nodes(j).value() > 0);
                              FEM_ASSERT(nodes(j).value() <= U.dim1());
                              Ue.xelem(3 * j + k) = U.xelem(nodes.xelem(j).value() - 1 + iNumNodesMesh * (k + 6 * l));
                         }
                    }

                    double dThetail = 0.;

                    if (dTheta.columns()) {
                         for (octave_idx_type j = 0; j < iNumNodes; ++j) {
                              dThetail += Ht.xelem(j) * dTheta.xelem(j + iNumNodes * l);
                         }
                    }

                    if (iNumPreLoads) {
                         const double epsilonTherm = gamma * dThetail;

                         for (octave_idx_type j = 0; j < iNumStrains; ++j) {
                              epsilonik.xelem(j) = j < 3 ? epsilonTherm : 0.;
                         }
                    }

                    if (epsilonRef.numel()) {
                         for (octave_idx_type j = 0; j < iNumNodes; ++j) {
                              for (octave_idx_type k = 0; k < iNumStrains; ++k) {
                                   epsilonik.xelem(k) += Ht.xelem(j) * epsilonRef.xelem(k + iNumStrains * (j + iNumNodes * l));
                              }
                         }
                    }

                    for (octave_idx_type j = 0; j < iNumStrains; ++j) {
                         double taugj = 0.;

                         for (octave_idx_type k = 0; k < iNumDof; ++k) {
                              taugj += CB.xelem(j + iNumStrains * k) * Ue.xelem(k);
                         }

                         if (iNumPreLoads) {
                              for (octave_idx_type k = 0; k < iNumStrains; ++k) {
                                   taugj -= C.xelem(j + iNumStrains * k) * epsilonik.xelem(k);
                              }
                         }

                         taug.xelem(i + iNumGauss * (l * iNumStrains + j)) = taugj;
                    }
               }
          }

          const Matrix tauen = InterpGaussToNodal(eMatType, taug);

          FEM_ASSERT(tauen.rows() == iNumNodes);
          FEM_ASSERT(tauen.columns() == taug.columns());

          for (octave_idx_type k = 0; k < iNumLoads; ++k) {
               for (octave_idx_type j = 0; j < iNumStrains; ++j) {
                    for (octave_idx_type i = 0; i < iNumNodes; ++i) {
                         taun.xelem(id - 1 + iNumElem * (i + iNumNodes * (j + k * iNumStrains))) = tauen.xelem(i + iNumNodes * (k * iNumStrains + j));
                    }
               }
          }
     }

     void ParticleVelocityGradient(const ColumnVector& rv,
                                   Matrix& vg,
                                   const octave_idx_type i,
                                   const octave_idx_type l,
                                   const Matrix& B,
                                   const ColumnVector& Phie,
                                   const ColumnVector& PhiPe,
                                   const double rho,
                                   const double tau,
                                   const octave_idx_type iNumDof) const {
          ParticleVelocityGradientNoPML<double>(rv, vg, i, l, B, Phie, PhiPe, rho, tau, iNumDof);
     }

     void ParticleVelocityGradient(const ColumnVector& rv,
                                   ComplexMatrix& vg,
                                   const octave_idx_type i,
                                   const octave_idx_type l,
                                   const Matrix& B,
                                   const ComplexColumnVector& Phie,
                                   const ComplexColumnVector& PhiPe,
                                   const double rho,
                                   const double tau,
                                   const octave_idx_type iNumDof) const {
          if (f.isempty()) {
               ParticleVelocityGradientNoPML<std::complex<double>>(rv, vg, i, l, B, Phie, PhiPe, rho, tau, iNumDof);
          } else {
               ParticleVelocityGradientPML(rv, vg, i, l, B, Phie, PhiPe, rho, tau, iNumDof);
          }
     }

     template <typename T>
     void ParticleVelocityGradientNoPML(const ColumnVector&,
                                        typename PostProcTypeTraits<T>::MatrixType& vg,
                                        const octave_idx_type i,
                                        const octave_idx_type l,
                                        const Matrix& B,
                                        const typename PostProcTypeTraits<T>::ColumnVectorType& Phie,
                                        const typename PostProcTypeTraits<T>::ColumnVectorType& PhiPe,
                                        const double rho,
                                        const double tau,
                                        const octave_idx_type iNumDof) const {
          constexpr octave_idx_type iNumComp = 3;
          const octave_idx_type iNumGauss = vg.rows();

          FEM_ASSERT(B.rows() == iNumComp);

          for (octave_idx_type j = 0; j < iNumComp; ++j) {
               T vgj{}, vgPj{};

               for (octave_idx_type k = 0; k < iNumDof; ++k) {
                    const double Bjk = B.xelem(j + iNumComp * k);
                    vgj += Bjk * Phie.xelem(k);
                    vgPj += Bjk * PhiPe.xelem(k);
               }

               vg.xelem(i + iNumGauss * (l * iNumComp + j)) = (vgj + tau * vgPj) / rho;
          }
     }

     void ParticleVelocityGradientPML(const ColumnVector& rv,
                                      ComplexMatrix& RFR_Tvg,
                                      const octave_idx_type i,
                                      const octave_idx_type l,
                                      const Matrix& B,
                                      const ComplexColumnVector& Phie,
                                      const ComplexColumnVector& PhiPe,
                                      const double rho,
                                      const double tau,
                                      const octave_idx_type iNumDof) const {
          constexpr octave_idx_type iNumComp = 3;
          const octave_idx_type iNumGauss = RFR_Tvg.rows();
          ComplexColumnVector vg(3);

          FEM_ASSERT(B.rows() == iNumComp);

          for (octave_idx_type j = 0; j < iNumComp; ++j) {
               std::complex<double> vgj{}, vgPj{};

               for (octave_idx_type k = 0; k < iNumDof; ++k) {
                    const double Bjk = B.xelem(j + iNumComp * k);
                    vgj += Bjk * Phie.xelem(k);
                    vgPj += Bjk * PhiPe.xelem(k);
               }

               vg.xelem(j) = (vgj + tau * vgPj) / rho;
          }

          Matrix R(3, 3), H(1, iNumDof);

          ScalarInterpMatrix(rv, H, 0);

          PMLCoordinateSystem(rv, H, R);

          ComplexColumnVector FR_Tvg(3);

          for (octave_idx_type j = 0; j < iNumComp; ++j) {
               std::complex<double> R_Tvgj{};

               for (octave_idx_type k = 0; k < iNumComp; ++k) {
                    R_Tvgj += R.xelem(k + 3 * j) * vg.xelem(k);
               }

               const auto fj = PMLF(j, H);

               FR_Tvg.xelem(j) = fj * R_Tvgj;
          }

          for (octave_idx_type j = 0; j < iNumComp; ++j) {
               std::complex<double> RFR_Tvgj{};

               for (octave_idx_type k = 0; k < 3; ++k) {
                    RFR_Tvgj += R.xelem(j + 3 * k) * FR_Tvg.xelem(k);
               }

               RFR_Tvg.xelem(i + iNumGauss * (l * iNumComp + j)) = RFR_Tvgj;
          }
     }

     template <typename T>
     void ParticleVelocityNodalElem(typename PostProcTypeTraits<T>::NDArrayType& vn,
                                    FemMatrixType eMatType,
                                    const typename PostProcTypeTraits<T>::NDArrayType& Phi,
                                    const typename PostProcTypeTraits<T>::NDArrayType& PhiP) const {
          typedef typename PostProcTypeTraits<T>::ColumnVectorType TColumnVector;
          typedef typename PostProcTypeTraits<T>::MatrixType TMatrix;

          const IntegrationRule& oIntegRule = GetIntegrationRule(eMatType);
          const octave_idx_type iNumGauss = oIntegRule.iGetNumEvalPoints();
          const octave_idx_type iNumDof = iGetNumDof(eMatType);
          const octave_idx_type iNumDir = oIntegRule.iGetNumDirections();
          const octave_idx_type iNumLoads = Phi.dim2();
          const octave_idx_type iNumNodes = nodes.numel();
          const octave_idx_type iNumNodesMesh = Phi.dim1();
          const octave_idx_type iNumElem = vn.dim1();
          constexpr octave_idx_type iNumComp = 3;

          ColumnVector rv(iNumDir);

          FEM_ASSERT(vn.ndims() >= 3);
          FEM_ASSERT(id >= 1);
          FEM_ASSERT(vn.dim1() >= id);
          FEM_ASSERT(vn.dim2() == iNumNodes);
          FEM_ASSERT(vn.dim3() == iNumComp);
          FEM_ASSERT((iNumLoads == 1 && vn.ndims() == 3) || (vn.dims()(3) == iNumLoads));
          FEM_ASSERT(Phi.ndims() >= 2);
          FEM_ASSERT(Phi.dim2() >= 1);
          FEM_ASSERT(Phi.dim1() == PhiP.dim1());
          FEM_ASSERT(iNumDof == X.columns());

          const double rho = material->Density();
          const double eta = material->ShearViscosity();
          const double zeta = material->VolumeViscosity();
          const double c = material->SpeedOfSound();
          const double tau = (4./3. * eta + zeta) / (rho * std::pow(c, 2));

          Matrix J(iNumDir, iNumDir), invJ(iNumDir, iNumDir), B(iNumComp, iNumDof);
          TColumnVector Phie(iNumDof), PhiPe(iNumDof);
          TMatrix vg(iNumGauss, iNumComp * iNumLoads);

          for (octave_idx_type i = 0; i < iNumGauss; ++i) {
               for (octave_idx_type j = 0; j < iNumDir; ++j) {
                    rv.xelem(j) = oIntegRule.dGetPosition(i, j);
               }

               const double detJ = Jacobian(rv, J);

               ScalarGradientMatrix(rv, J, detJ, invJ, B);

               for (octave_idx_type l = 0; l < iNumLoads; ++l) {
                    for (octave_idx_type j = 0; j < iNumNodes; ++j) {
                         const octave_idx_type idof = nodes.xelem(j).value() - 1;

                         Phie.xelem(j) = Phi.xelem(idof + iNumNodesMesh * l);
                         PhiPe.xelem(j) = PhiP.xelem(idof + iNumNodesMesh * l);
                    }

                    ParticleVelocityGradient(rv, vg, i, l, B, Phie, PhiPe, rho, tau, iNumDof);
               }
          }

          const TMatrix ven = InterpGaussToNodal(eMatType, vg);

          FEM_ASSERT(ven.rows() == iNumNodes);
          FEM_ASSERT(ven.columns() == vg.columns());

          for (octave_idx_type k = 0; k < iNumLoads; ++k) {
               for (octave_idx_type j = 0; j < iNumComp; ++j) {
                    for (octave_idx_type i = 0; i < iNumNodes; ++i) {
                         vn.xelem(id - 1 + iNumElem * (i + iNumNodes * (j + k * iNumComp))) = ven.xelem(i + iNumNodes * (k * iNumComp + j));
                    }
               }
          }
     }

     void ThermalLoadVector(Matrix& R, MeshInfo& info, FemMatrixType eMatType) const {
          const IntegrationRule& oIntegRule = GetIntegrationRule(eMatType);
          const octave_idx_type iNumGauss = oIntegRule.iGetNumEvalPoints();
          const octave_idx_type iNumDir = oIntegRule.iGetNumDirections();
          const octave_idx_type iNumNodes = nodes.numel();
          const octave_idx_type iNumDof = iGetNumDof(eMatType);
          const Matrix& C = material->LinearElasticity();
          const double gamma = material->ThermalExpansion();
          const octave_idx_type iNumStrains = C.rows();

          FEM_ASSERT(R.rows() == iNumDof);
          FEM_ASSERT(R.columns() == iNumPreLoads);

          Matrix Ht(1, iNumNodes);
          ColumnVector rv(iNumDir);

          Matrix J(iNumDir, iNumDir), invJ(iNumDir, iNumDir), B(iNumStrains, iNumDof);
          Matrix epsilon(iNumStrains, iNumPreLoads);
          Matrix Cepsilon(iNumStrains, iNumPreLoads);

          for (octave_idx_type i = 0; i < iNumGauss; ++i) {
               const double alpha = oIntegRule.dGetWeight(i);

               for (octave_idx_type j = 0; j < iNumDir; ++j) {
                    rv.xelem(j) = oIntegRule.dGetPosition(i, j);
               }

               ScalarInterpMatrix(rv, Ht, 0);

               for (octave_idx_type k = 0; k < iNumPreLoads; ++k) {
                    double dThetaik = 0.;

                    if (dTheta.columns()) {
                         for (octave_idx_type j = 0; j < iNumNodes; ++j) {
                              dThetaik += Ht.xelem(j) * dTheta.xelem(j + iNumNodes * k);
                         }
                    }

                    FEM_ASSERT(iNumStrains >= 3);

                    const double epsilonk = gamma * dThetaik;

                    for (octave_idx_type j = 0; j < iNumStrains; ++j) {
                         epsilon.xelem(j + iNumStrains * k) = j < 3 ? epsilonk : 0.;
                    }

                    if (epsilonRef.numel()) {
                         for (octave_idx_type j = 0; j < iNumNodes; ++j) {
                              for (octave_idx_type l = 0; l < iNumStrains; ++l) {
                                   epsilon.xelem(l + iNumStrains * k) += Ht.xelem(j) * epsilonRef.xelem(l + iNumStrains * (j + iNumNodes * k));
                              }
                         }
                    }
               }

               for (octave_idx_type l = 0; l < iNumPreLoads; ++l) {
                    for (octave_idx_type j = 0; j < iNumStrains; ++j) {
                         double Cepsilonjl = 0.;

                         for (octave_idx_type k = 0; k < iNumStrains; ++k) {
                              Cepsilonjl += C.xelem(j + iNumStrains * k) * epsilon.xelem(k + iNumStrains * l);
                         }

                         Cepsilon.xelem(j + iNumStrains * l) = Cepsilonjl;
                    }
               }

               const double detJ = Jacobian(rv, J);

               AddMeshInfo(info, oIntegRule, detJ);

               StrainMatrix(rv, J, detJ, invJ, B);

               for (octave_idx_type j = 0; j < iNumPreLoads; ++j) {
                    for (octave_idx_type k = 0; k < iNumDof; ++k) {
                         double Rkj = 0.;

                         for (octave_idx_type l = 0; l < iNumStrains; ++l) {
                              Rkj += B.xelem(l + iNumStrains * k) * Cepsilon.xelem(l + iNumStrains * j);
                         }

                         R.xelem(k + iNumDof * j) += detJ * alpha * Rkj;
                    }
               }
          }
     }

     void AcousticStiffnessMatrixCSL(const ColumnVector&, Matrix& k, FemMatrixType eMatType) const {
          const double sign = (eMatType == MAT_STIFFNESS_ACOUSTICS_RE ? 1. : -1.);
          const double coef = sign / material->Density();

          FEM_ASSERT(k.rows() == 3);
          FEM_ASSERT(k.columns() == 3);

          for (octave_idx_type i = 0; i < 3; ++i) {
               k.xelem(i + 3 * i) = coef;
          }
     }

     void AcousticStiffnessMatrix(Matrix& Ke, MeshInfo& info, FemMatrixType eMatType) const {
          ScalarFieldStiffnessMatrix(&Element3D::AcousticStiffnessMatrixCSL, Ke, info, eMatType);
     }

     void PMLCoordinateSystem(const ColumnVector& rv, const Matrix& H, Matrix& R) const {
          FEM_ASSERT(R.rows() == 3);
          FEM_ASSERT(R.columns() == 3);
          FEM_ASSERT(H.rows() == 1);
          FEM_ASSERT(H.columns() == nodes.numel());
          FEM_ASSERT(e1.rows() == 3);
          FEM_ASSERT(e2.rows() == 3);
          FEM_ASSERT(e1.columns() == nodes.numel());
          FEM_ASSERT(e2.columns() == nodes.numel());

          for (octave_idx_type i = 0; i < 3; ++i) {
               double e1i = 0., e2i = 0.;

               for (octave_idx_type k = 0; k < H.columns(); ++k) {
                    const double Hk = H.xelem(k);

                    e1i += Hk * e1.xelem(i + 3 * k);
                    e2i += Hk * e2.xelem(i + 3 * k);
               }

               R.xelem(i + 3 * 0) = e1i;
               R.xelem(i + 3 * 1) = e2i;
          }

          R.xelem(0 + 3 * 2) = R.xelem(1 + 3 * 0) * R.xelem(2 + 3 * 1) - R.xelem(2 + 3 * 0) * R.xelem(1 + 3 * 1);
          R.xelem(1 + 3 * 2) = R.xelem(2 + 3 * 0) * R.xelem(0 + 3 * 1) - R.xelem(0 + 3 * 0) * R.xelem(2 + 3 * 1);
          R.xelem(2 + 3 * 2) = R.xelem(0 + 3 * 0) * R.xelem(1 + 3 * 1) - R.xelem(1 + 3 * 0) * R.xelem(0 + 3 * 1);

          R.xelem(0 + 3 * 1) = R.xelem(1 + 3 * 2) * R.xelem(2 + 3 * 0) - R.xelem(2 + 3 * 2) * R.xelem(1 + 3 * 0);
          R.xelem(1 + 3 * 1) = R.xelem(2 + 3 * 2) * R.xelem(0 + 3 * 0) - R.xelem(0 + 3 * 2) * R.xelem(2 + 3 * 0);
          R.xelem(2 + 3 * 1) = R.xelem(0 + 3 * 2) * R.xelem(1 + 3 * 0) - R.xelem(1 + 3 * 2) * R.xelem(0 + 3 * 0);

          for (octave_idx_type j = 0; j < 3; ++j) {
               double nej = 0.;

               for (octave_idx_type i = 0; i < 3; ++i) {
                    const double Rij = R.xelem(i + 3 * j);
                    nej += Rij * Rij;
               }

               if (nej <= 0.) {
                    const ElementTypes::TypeInfo& oInfo = ElementTypes::GetType(eltype);
                    throw std::runtime_error("element 3D: vectors mesh.elements.perfectly_matched_layers."s
                                             + oInfo.name
                                             + ".e1 and mesh.elements.perfectly_matched_layers."
                                             + oInfo.name
                                             + ".e2 are co-linear");
               }

               nej = sqrt(nej);

               for (octave_idx_type i = 0; i < 3; ++i) {
                    R.xelem(i + 3 * j) /= nej;
               }
          }
     }

     std::complex<double> PMLF(const octave_idx_type i, const Matrix& H) const {
          std::complex<double> fi{};

          FEM_ASSERT(f.rows() == 3);
          FEM_ASSERT(f.columns() == nodes.numel());

          for (octave_idx_type j = 0; j < H.columns(); ++j) {
               fi += f.xelem(i + 3 * j) * H.xelem(j);
          }

          return fi;
     }

     void AcousticStiffnessDampingMatrixPMLCSL(const ColumnVector& rv, Matrix& RF_2R_T, FemMatrixType eMatType) const {
          // Apply transformation for perfectly matched layers

          const octave_idx_type iNumNodes = nodes.numel();

          FEM_ASSERT(RF_2R_T.rows() == 3);
          FEM_ASSERT(RF_2R_T.columns() == 3);

          Matrix H(1, iNumNodes);

          ScalarInterpMatrix(rv, H, 0);

          Matrix R(3, 3);

          PMLCoordinateSystem(rv, H, R);

          double sign;

          switch (eMatType) {
          case MAT_STIFFNESS_ACOUSTICS_RE:
          case MAT_STIFFNESS_ACOUSTICS_IM:
          case MAT_DAMPING_ACOUSTICS_RE:
          case MAT_DAMPING_ACOUSTICS_IM:
               sign = 1.;
               break;
          case MAT_STIFFNESS_FLUID_STRUCT_RE:
          case MAT_STIFFNESS_FLUID_STRUCT_IM:
          case MAT_DAMPING_FLUID_STRUCT_RE:
          case MAT_DAMPING_FLUID_STRUCT_IM:
               sign = -1.;
               break;
          default:
               FEM_ASSERT(false);
               throw std::logic_error("element 3D: invalid matrix type");
          }

          double coef;

          switch (eMatType) {
          case MAT_STIFFNESS_ACOUSTICS_RE:
          case MAT_STIFFNESS_ACOUSTICS_IM:
          case MAT_STIFFNESS_FLUID_STRUCT_RE:
          case MAT_STIFFNESS_FLUID_STRUCT_IM:
               coef = sign / material->Density();
               break;
          case MAT_DAMPING_ACOUSTICS_RE:
          case MAT_DAMPING_ACOUSTICS_IM:
          case MAT_DAMPING_FLUID_STRUCT_RE:
          case MAT_DAMPING_FLUID_STRUCT_IM: {
               const double rho = material->Density();
               const double eta = material->ShearViscosity();
               const double zeta = material->VolumeViscosity();
               const double c = material->SpeedOfSound();
               coef = sign * (4./3. * eta + zeta) / std::pow(rho * c, 2);
          } break;
          default:
               FEM_ASSERT(false);
               throw std::logic_error("element 3D: invalid matrix type");
          }

          ComplexColumnVector frv(3);

          std::complex<double> detJ{1., 0.};

          for (octave_idx_type i = 0; i < 3; ++i) {
               detJ /= (frv.xelem(i) = PMLF(i, H));
          }

          ColumnVector F_2(3);

          for (octave_idx_type i = 0; i < 3; ++i) {
               double fi_2_detJ;

               switch (eMatType) {
               case MAT_STIFFNESS_ACOUSTICS_RE:
               case MAT_DAMPING_ACOUSTICS_RE:
               case MAT_STIFFNESS_FLUID_STRUCT_RE:
               case MAT_DAMPING_FLUID_STRUCT_RE:
                    //realpart(f^2*detJ);
                    //detJRe*(fxRe^2-fxIm^2)-2*detJIm*fxIm*fxRe
                    fi_2_detJ = std::real(detJ) * (std::pow(std::real(frv.xelem(i)), 2) - std::pow(std::imag(frv.xelem(i)), 2))
                              - 2. * std::imag(detJ) * std::imag(frv.xelem(i)) * std::real(frv.xelem(i));
                    break;
               case MAT_STIFFNESS_ACOUSTICS_IM:
               case MAT_DAMPING_ACOUSTICS_IM:
               case MAT_STIFFNESS_FLUID_STRUCT_IM:
               case MAT_DAMPING_FLUID_STRUCT_IM:
                    //imagpart(fx^2*detJ);
                    //detJIm*(fxRe^2-fxIm^2)+2*detJRe*fxIm*fxRe
                    fi_2_detJ = std::imag(detJ) * (std::pow(std::real(frv.xelem(i)), 2) - std::pow(std::imag(frv.xelem(i)), 2))
                              + 2. * std::real(detJ) * std::imag(frv.xelem(i)) * std::real(frv.xelem(i));
                    break;
               default:
                    FEM_ASSERT(false);
                    throw std::logic_error("element 3D: invalid matrix type");
               }

               F_2.xelem(i) = coef * fi_2_detJ;
          }

          Matrix F_2R_T(3, 3);

          for (octave_idx_type j = 0; j < 3; ++j) {
               for (octave_idx_type i = 0; i < 3; ++i) {
                    F_2R_T.xelem(i + 3 * j) = R.xelem(j + 3 * i) * F_2.xelem(i);
               }
          }

          for (octave_idx_type j = 0; j < 3; ++j) {
               for (octave_idx_type i = 0; i < 3; ++i) {
                    double RF_2R_Tij = 0.;

                    for (octave_idx_type k = 0; k < 3; ++k) {
                         RF_2R_Tij += R.xelem(i + 3 * k) * F_2R_T.xelem(k + 3 * j);
                    }

                    RF_2R_T.xelem(i + 3 * j) = RF_2R_Tij;
               }
          }
     }

     void AcousticStiffnessMatrixPML(Matrix& Ke, MeshInfo& info, FemMatrixType eMatType) const {
          ScalarFieldStiffnessMatrix(&Element3D::AcousticStiffnessDampingMatrixPMLCSL, Ke, info, eMatType);
     }

     void AcousticDampingMatrixPML(Matrix& De, MeshInfo& info, FemMatrixType eMatType) const {
          ScalarFieldStiffnessMatrix(&Element3D::AcousticStiffnessDampingMatrixPMLCSL, De, info, eMatType);
     }

     void AcousticDampingMatrixCSL(const ColumnVector&, Matrix& k, FemMatrixType eMatType) const {
          const double eta = material->ShearViscosity();
          const double zeta = material->VolumeViscosity();
          const double rho = material->Density();
          const double c = material->SpeedOfSound();
          const double sign = (eMatType == MAT_DAMPING_ACOUSTICS_RE ? 1. : -1.);
          const double coef = sign * (4./3. * eta + zeta) / std::pow(rho * c, 2);

          FEM_ASSERT(k.rows() == 3);
          FEM_ASSERT(k.columns() == 3);

          for (octave_idx_type i = 0; i < 3; ++i) {
               k.xelem(i + 3 * i) = coef;
          }
     }

     void AcousticDampingMatrix(Matrix& De, MeshInfo& info, FemMatrixType eMatType) const {
          ScalarFieldStiffnessMatrix(&Element3D::AcousticDampingMatrixCSL, De, info, eMatType);
     }

     void ThermalConductivityMatrixCSL(const ColumnVector&, Matrix& k, FemMatrixType) const {
          FEM_ASSERT(material->ThermalConductivity().rows() == 3);
          FEM_ASSERT(material->ThermalConductivity().columns() == 3);

          k = material->ThermalConductivity();
     }

     void ThermalConductivityMatrix(Matrix& Ke, MeshInfo& info, FemMatrixType eMatType) const {
          ScalarFieldStiffnessMatrix(&Element3D::ThermalConductivityMatrixCSL, Ke, info, eMatType);
     }

     typedef void (Element3D::*ScalarFieldCSL)(const ColumnVector&, Matrix&, FemMatrixType) const;

     void ScalarFieldStiffnessMatrix(ScalarFieldCSL pfnk, Matrix& Ke, MeshInfo& info, FemMatrixType eMatType) const {
          const IntegrationRule& oIntegRule = GetIntegrationRule(eMatType);
          const octave_idx_type iNumDof = iGetNumDof(eMatType);
          const octave_idx_type iNumDir = oIntegRule.iGetNumDirections();
          ColumnVector rv(iNumDir);
          constexpr octave_idx_type iNumStrains = 3;

          Matrix k(iNumStrains, iNumStrains, 0.);

          Matrix J(iNumDir, iNumDir), invJ(iNumDir, iNumDir), B(iNumStrains, iNumDof), kB(iNumStrains, iNumDof);

          for (octave_idx_type i = 0; i < oIntegRule.iGetNumEvalPoints(); ++i) {
               const double alpha = oIntegRule.dGetWeight(i);

               for (octave_idx_type j = 0; j < iNumDir; ++j) {
                    rv.xelem(j) = oIntegRule.dGetPosition(i, j);
               }

               const double detJ = Jacobian(rv, J);

               AddMeshInfo(info, oIntegRule, detJ);

               ScalarGradientMatrix(rv, J, detJ, invJ, B);

               (this->*pfnk)(rv, k, eMatType);

               for (octave_idx_type l = 0; l < iNumStrains; ++l) {
                    for (octave_idx_type m = 0; m < iNumDof; ++m) {
                         double kBlm = 0.;

                         for (octave_idx_type n = 0; n < iNumStrains; ++n) {
                              kBlm += detJ * alpha * k.xelem(l + 3 * n) * B.xelem(n + iNumStrains * m);
                         }

                         FEM_ASSERT(std::isfinite(kBlm));

                         kB.xelem(l + iNumStrains * m) = kBlm;
                    }
               }

#ifdef DEBUG
               for (octave_idx_type k = 0; k < kB.rows(); ++k) {
                    for (octave_idx_type j = 0; j < kB.columns(); ++j) {
                         FEM_ASSERT(std::isfinite(kB(k, j)));
                    }
               }
#endif

               for (octave_idx_type l = 0; l < iNumDof; ++l) {
                    for (octave_idx_type m = l; m < iNumDof; ++m) {
                         double Kelm = 0.;

                         for (octave_idx_type n = 0; n < iNumStrains; ++n) {
                              Kelm += B.xelem(n + iNumStrains * l) * kB.xelem(n + iNumStrains * m);
                         }

                         FEM_ASSERT(std::isfinite(Kelm));

                         Ke.xelem(l + iNumDof * m) += Kelm;
                    }
               }
          }

          for (octave_idx_type i = 1; i < iNumDof; ++i) {
               for (octave_idx_type j = 0; j < i; ++j) {
                    Ke.xelem(i + iNumDof * j) = Ke.xelem(j + iNumDof * i);
               }
          }

#ifdef DEBUG
          for (octave_idx_type i = 0; i < Ke.rows(); ++i) {
               for (octave_idx_type j = 0; j < Ke.columns(); ++j) {
                    FEM_ASSERT(std::isfinite(Ke(i, j)));
               }
          }
#endif
     }

     typedef double (Element3D::*ScalarFieldMassCSL)(const ColumnVector&, const Matrix&, FemMatrixType) const;

     double AcousticMassCSL(const ColumnVector&, const Matrix&, FemMatrixType eMatType) const {
          double sign;

          switch (eMatType) {
          case MAT_MASS_ACOUSTICS_RE:
          case MAT_MASS_ACOUSTICS_IM:
               sign = 1.;
               break;
          case MAT_MASS_FLUID_STRUCT_RE:
          case MAT_MASS_FLUID_STRUCT_IM:
               sign = -1.;
               break;
          default:
               throw std::logic_error("element 3D: invalid matrix type");
          }

          return sign / (material->Density() * std::pow(material->SpeedOfSound(), 2));
     }

     double AcousticMassCSLPML(const ColumnVector& rv, const Matrix& H, FemMatrixType eMatType) const {
          std::complex<double> detJ{1., 0.};

          for (octave_idx_type i = 0; i < 3; ++i) {
               detJ /= PMLF(i, H);
          }

          const double coef = AcousticMassCSL(rv, H, eMatType);

          switch (eMatType) {
          case MAT_MASS_ACOUSTICS_RE:
          case MAT_MASS_FLUID_STRUCT_RE:
               return std::real(detJ) * coef;
          case MAT_MASS_ACOUSTICS_IM:
          case MAT_MASS_FLUID_STRUCT_IM:
               return std::imag(detJ) * coef;
          default:
               throw std::logic_error("element 3D: invalid matrix type");
          }
     }

     void AcousticMassMatrix(Matrix& Ke, MeshInfo& info, FemMatrixType eMatType) const {
          ScalarFieldMassMatrix(&Element3D::AcousticMassCSL, Ke, info, eMatType);
     }

     void AcousticMassMatrixPML(Matrix& Ke, MeshInfo& info, FemMatrixType eMatType) const {
          ScalarFieldMassMatrix(&Element3D::AcousticMassCSLPML, Ke, info, eMatType);
     }

     double HeatCapacityCSL(const ColumnVector&, const Matrix&, FemMatrixType) const {
          return material->Density() * material->HeatCapacity();
     }

     void HeatCapacityMatrix(Matrix& Ce, MeshInfo& info, FemMatrixType eMatType) const {
          ScalarFieldMassMatrix(&Element3D::HeatCapacityCSL, Ce, info, eMatType);
     }

     void ScalarFieldMassMatrix(ScalarFieldMassCSL pfncoef, Matrix& Ce, MeshInfo& info, FemMatrixType eMatType) const {
          const IntegrationRule& oIntegRule = GetIntegrationRule(eMatType);
          const octave_idx_type iNumDof = iGetNumDof(eMatType);
          const octave_idx_type iNumDir = oIntegRule.iGetNumDirections();
          ColumnVector rv(iNumDir);
          Matrix J(iNumDir, iNumDir), H(1, iNumDof);

          for (octave_idx_type i = 0; i < oIntegRule.iGetNumEvalPoints(); ++i) {
               const double alpha = oIntegRule.dGetWeight(i);

               for (octave_idx_type j = 0; j < iNumDir; ++j) {
                    rv.xelem(j) = oIntegRule.dGetPosition(i, j);
               }

               const double detJ = Jacobian(rv, J);

               AddMeshInfo(info, oIntegRule, detJ);

               ScalarInterpMatrix(rv, H, 0);

               const double coef = (this->*pfncoef)(rv, H, eMatType);

               for (octave_idx_type l = 0; l < iNumDof; ++l) {
                    for (octave_idx_type m = l; m < iNumDof; ++m) {
                         Ce.xelem(l + iNumDof * m) += H.xelem(l) * H.xelem(m) * alpha * coef * detJ;
                    }
               }
          }

          for (octave_idx_type i = 1; i < iNumDof; ++i) {
               for (octave_idx_type j = 0; j < i; ++j) {
                    Ce.xelem(i + iNumDof * j) = Ce.xelem(j + iNumDof * i);
               }
          }
     }

     void CollocationPoints(NDArray& Xg, FemMatrixType eMatType) const {
          const IntegrationRule& oIntegRule = GetIntegrationRule(eMatType);
          const octave_idx_type iNumNodes = nodes.numel();
          const octave_idx_type iNumDir = oIntegRule.iGetNumDirections();
          const octave_idx_type iNumDisp = X.rows();
          const octave_idx_type iNumElem = Xg.dim1();
          const octave_idx_type iNumGauss = oIntegRule.iGetNumEvalPoints();

          FEM_ASSERT(Xg.ndims() == 3);
          FEM_ASSERT(Xg.rows() >= id);
          FEM_ASSERT(Xg.columns() == oIntegRule.iGetNumEvalPoints());
          FEM_ASSERT(Xg.pages() == iNumDisp);

          Matrix H(1, iNumNodes);
          ColumnVector rv(iNumDir);

          for (octave_idx_type i = 0; i < iNumGauss; ++i) {
               for (octave_idx_type j = 0; j < iNumDir; ++j) {
                    rv.xelem(j) = oIntegRule.dGetPosition(i, j);
               }

               ScalarInterpMatrix(rv, H, 0);

               for (octave_idx_type k = 0; k < iNumDisp; ++k) {
                    double Xk = 0.;

                    for (octave_idx_type l = 0; l < iNumNodes; ++l) {
                         Xk += H.xelem(l) * X.xelem(k + iNumDisp * l);
                    }

                    Xg.xelem(id - 1 + iNumElem * (i + iNumGauss * k)) = Xk;
               }
          }
     }

     template <typename T>
     static T Determinant3x3(const Array<T>& J) {
          FEM_ASSERT(J.rows() == 3);
          FEM_ASSERT(J.columns() == 3);

          return J.xelem(0)*(J.xelem(4)*J.xelem(8)-J.xelem(5)*J.xelem(7))-J.xelem(3)*(J.xelem(1)*J.xelem(8)-J.xelem(2)*J.xelem(7))+(J.xelem(1)*J.xelem(5)-J.xelem(2)*J.xelem(4))*J.xelem(6);
     }

     template <typename T>
     static void Inverse3x3(const Array<T>& J, const T detJ, Array<T>& invJ) {
          FEM_ASSERT(J.rows() == 3);
          FEM_ASSERT(J.columns() == 3);
          FEM_ASSERT(invJ.rows() == 3);
          FEM_ASSERT(invJ.columns() == 3);

          invJ.xelem(0) = (J.xelem(4)*J.xelem(8)-J.xelem(5)*J.xelem(7))/detJ;
          invJ.xelem(1) = -(J.xelem(1)*J.xelem(8)-J.xelem(2)*J.xelem(7))/detJ;
          invJ.xelem(2) = (J.xelem(1)*J.xelem(5)-J.xelem(2)*J.xelem(4))/detJ;
          invJ.xelem(3) = -(J.xelem(3)*J.xelem(8)-J.xelem(5)*J.xelem(6))/detJ;
          invJ.xelem(4) = (J.xelem(0)*J.xelem(8)-J.xelem(2)*J.xelem(6))/detJ;
          invJ.xelem(5) = -(J.xelem(0)*J.xelem(5)-J.xelem(2)*J.xelem(3))/detJ;
          invJ.xelem(6) = (J.xelem(3)*J.xelem(7)-J.xelem(4)*J.xelem(6))/detJ;
          invJ.xelem(7) = -(J.xelem(0)*J.xelem(7)-J.xelem(1)*J.xelem(6))/detJ;
          invJ.xelem(8) = (J.xelem(0)*J.xelem(4)-J.xelem(1)*J.xelem(3))/detJ;
     }

     template <typename T>
     static T Determinant4x4(const Array<T>& J) {
          FEM_ASSERT(J.rows() == 4);
          FEM_ASSERT(J.columns() == 4);

          return J.xelem(0)*(J.xelem(5)*(J.xelem(10)*J.xelem(15)-J.xelem(11)*J.xelem(14))-J.xelem(9)*(J.xelem(6)*J.xelem(15)-J.xelem(7)*J.xelem(14))+(J.xelem(6)*J.xelem(11)-J.xelem(7)*J.xelem(10))*J.xelem(13))-J.xelem(4)*(J.xelem(1)*(J.xelem(10)*J.xelem(15)-J.xelem(11)*J.xelem(14))-J.xelem(9)*(J.xelem(2)*J.xelem(15)-J.xelem(3)*J.xelem(14))+(J.xelem(2)*J.xelem(11)-J.xelem(3)*J.xelem(10))*J.xelem(13))+J.xelem(8)*(J.xelem(1)*(J.xelem(6)*J.xelem(15)-J.xelem(7)*J.xelem(14))-J.xelem(5)*(J.xelem(2)*J.xelem(15)-J.xelem(3)*J.xelem(14))+(J.xelem(2)*J.xelem(7)-J.xelem(3)*J.xelem(6))*J.xelem(13))-(J.xelem(1)*(J.xelem(6)*J.xelem(11)-J.xelem(7)*J.xelem(10))-J.xelem(5)*(J.xelem(2)*J.xelem(11)-J.xelem(3)*J.xelem(10))+(J.xelem(2)*J.xelem(7)-J.xelem(3)*J.xelem(6))*J.xelem(9))*J.xelem(12);
     }

     template <typename T>
     static void Inverse4x4(const Array<T>& J, const T detJ, Array<T>& invJ) {
          FEM_ASSERT(J.rows() == 4);
          FEM_ASSERT(J.columns() == 4);
          FEM_ASSERT(invJ.rows() == 4);
          FEM_ASSERT(invJ.columns() == 4);

          invJ.xelem(0) = ((J.xelem(5)*J.xelem(10)-J.xelem(6)*J.xelem(9))*J.xelem(15)+(J.xelem(7)*J.xelem(9)-J.xelem(5)*J.xelem(11))*J.xelem(14)+(J.xelem(6)*J.xelem(11)-J.xelem(7)*J.xelem(10))*J.xelem(13))/detJ;
          invJ.xelem(1) = -((J.xelem(1)*J.xelem(10)-J.xelem(2)*J.xelem(9))*J.xelem(15)+(J.xelem(3)*J.xelem(9)-J.xelem(1)*J.xelem(11))*J.xelem(14)+(J.xelem(2)*J.xelem(11)-J.xelem(3)*J.xelem(10))*J.xelem(13))/detJ;
          invJ.xelem(2) = ((J.xelem(1)*J.xelem(6)-J.xelem(2)*J.xelem(5))*J.xelem(15)+(J.xelem(3)*J.xelem(5)-J.xelem(1)*J.xelem(7))*J.xelem(14)+(J.xelem(2)*J.xelem(7)-J.xelem(3)*J.xelem(6))*J.xelem(13))/detJ;
          invJ.xelem(3) = -((J.xelem(1)*J.xelem(6)-J.xelem(2)*J.xelem(5))*J.xelem(11)+(J.xelem(3)*J.xelem(5)-J.xelem(1)*J.xelem(7))*J.xelem(10)+(J.xelem(2)*J.xelem(7)-J.xelem(3)*J.xelem(6))*J.xelem(9))/detJ;
          invJ.xelem(4) = -((J.xelem(4)*J.xelem(10)-J.xelem(6)*J.xelem(8))*J.xelem(15)+(J.xelem(7)*J.xelem(8)-J.xelem(4)*J.xelem(11))*J.xelem(14)+(J.xelem(6)*J.xelem(11)-J.xelem(7)*J.xelem(10))*J.xelem(12))/detJ;
          invJ.xelem(5) = ((J.xelem(0)*J.xelem(10)-J.xelem(2)*J.xelem(8))*J.xelem(15)+(J.xelem(3)*J.xelem(8)-J.xelem(0)*J.xelem(11))*J.xelem(14)+(J.xelem(2)*J.xelem(11)-J.xelem(3)*J.xelem(10))*J.xelem(12))/detJ;
          invJ.xelem(6) = -((J.xelem(0)*J.xelem(6)-J.xelem(2)*J.xelem(4))*J.xelem(15)+(J.xelem(3)*J.xelem(4)-J.xelem(0)*J.xelem(7))*J.xelem(14)+(J.xelem(2)*J.xelem(7)-J.xelem(3)*J.xelem(6))*J.xelem(12))/detJ;
          invJ.xelem(7) = ((J.xelem(0)*J.xelem(6)-J.xelem(2)*J.xelem(4))*J.xelem(11)+(J.xelem(3)*J.xelem(4)-J.xelem(0)*J.xelem(7))*J.xelem(10)+(J.xelem(2)*J.xelem(7)-J.xelem(3)*J.xelem(6))*J.xelem(8))/detJ;
          invJ.xelem(8) = ((J.xelem(4)*J.xelem(9)-J.xelem(5)*J.xelem(8))*J.xelem(15)+(J.xelem(7)*J.xelem(8)-J.xelem(4)*J.xelem(11))*J.xelem(13)+(J.xelem(5)*J.xelem(11)-J.xelem(7)*J.xelem(9))*J.xelem(12))/detJ;
          invJ.xelem(9) = -((J.xelem(0)*J.xelem(9)-J.xelem(1)*J.xelem(8))*J.xelem(15)+(J.xelem(3)*J.xelem(8)-J.xelem(0)*J.xelem(11))*J.xelem(13)+(J.xelem(1)*J.xelem(11)-J.xelem(3)*J.xelem(9))*J.xelem(12))/detJ;
          invJ.xelem(10) = ((J.xelem(0)*J.xelem(5)-J.xelem(1)*J.xelem(4))*J.xelem(15)+(J.xelem(3)*J.xelem(4)-J.xelem(0)*J.xelem(7))*J.xelem(13)+(J.xelem(1)*J.xelem(7)-J.xelem(3)*J.xelem(5))*J.xelem(12))/detJ;
          invJ.xelem(11) = -((J.xelem(0)*J.xelem(5)-J.xelem(1)*J.xelem(4))*J.xelem(11)+(J.xelem(3)*J.xelem(4)-J.xelem(0)*J.xelem(7))*J.xelem(9)+(J.xelem(1)*J.xelem(7)-J.xelem(3)*J.xelem(5))*J.xelem(8))/detJ;
          invJ.xelem(12) = -((J.xelem(4)*J.xelem(9)-J.xelem(5)*J.xelem(8))*J.xelem(14)+(J.xelem(6)*J.xelem(8)-J.xelem(4)*J.xelem(10))*J.xelem(13)+(J.xelem(5)*J.xelem(10)-J.xelem(6)*J.xelem(9))*J.xelem(12))/detJ;
          invJ.xelem(13) = ((J.xelem(0)*J.xelem(9)-J.xelem(1)*J.xelem(8))*J.xelem(14)+(J.xelem(2)*J.xelem(8)-J.xelem(0)*J.xelem(10))*J.xelem(13)+(J.xelem(1)*J.xelem(10)-J.xelem(2)*J.xelem(9))*J.xelem(12))/detJ;
          invJ.xelem(14) = -((J.xelem(0)*J.xelem(5)-J.xelem(1)*J.xelem(4))*J.xelem(14)+(J.xelem(2)*J.xelem(4)-J.xelem(0)*J.xelem(6))*J.xelem(13)+(J.xelem(1)*J.xelem(6)-J.xelem(2)*J.xelem(5))*J.xelem(12))/detJ;
          invJ.xelem(15) = ((J.xelem(0)*J.xelem(5)-J.xelem(1)*J.xelem(4))*J.xelem(10)+(J.xelem(2)*J.xelem(4)-J.xelem(0)*J.xelem(6))*J.xelem(9)+(J.xelem(1)*J.xelem(6)-J.xelem(2)*J.xelem(5))*J.xelem(8))/detJ;
     }

     virtual double Jacobian(const ColumnVector& rv, Matrix& J) const {
          FEM_ASSERT(J.rows() == 3);
          FEM_ASSERT(J.columns() == 3);
          FEM_ASSERT(rv.numel() == 3);

          const octave_idx_type N = X.columns();

          FEM_ASSERT(N == nodes.numel());
          FEM_ASSERT(X.rows() == 3);

          Matrix Hd(N, 3);

          ScalarInterpMatrixDer(rv, Hd);

          for (octave_idx_type i = 0; i < 3; ++i) {
               for (octave_idx_type j = 0; j < 3; ++j) {
                    double Jji = 0.;

                    for (octave_idx_type k = 0; k < N; ++k) {
                         Jji += X.xelem(i + 3 * k) * Hd.xelem(k + j * N);
                    }

                    J.xelem(j + i * 3) = Jji;
               }
          }

          return Determinant3x3(J);
     }

     virtual void StrainMatrix(const ColumnVector& rv, const Matrix& J, double detJ, Matrix& invJ, Matrix& B) const {
          FEM_ASSERT(rv.numel() == 3);
          FEM_ASSERT(J.rows() == 3);
          FEM_ASSERT(J.columns() == 3);
          FEM_ASSERT(invJ.rows() == 3);
          FEM_ASSERT(invJ.columns() == 3);
          FEM_ASSERT(B.rows() == 6);
          FEM_ASSERT(B.columns() == X.columns() * 3);

          Inverse3x3(J, detJ, invJ);

          const octave_idx_type N = X.columns();

          Matrix Hd(N, 3);

          ScalarInterpMatrixDer(rv, Hd);

          Matrix H0d(N, 3);

          for (octave_idx_type i = 0; i < N; ++i) {
               for (octave_idx_type j = 0; j < 3; ++j) {
                    double H0dij = 0.;

                    for (octave_idx_type k = 0; k < 3; ++k) {
                         H0dij += Hd.xelem(i + k * N) * invJ.xelem(j + k * 3);
                    }

                    H0d.xelem(i + j * N) = H0dij;
               }
          }

          static constexpr octave_idx_type i1[] = {0, 1, 2, 3, 3, 4, 4, 5, 5};
          static constexpr octave_idx_type i2[] = {0, 1, 2, 0, 1, 1, 2, 0, 2};
          static constexpr octave_idx_type i3[] = {0, 1, 2, 1, 0, 2, 1, 2, 0};

          static constexpr octave_idx_type M = sizeof(i1) / sizeof(i1[0]);

          for (octave_idx_type i = 0; i < 6 * N; ++i) {
               B.xelem(i) = 0.;
          }

          for (octave_idx_type k = 0; k < N; ++k) {
               for (octave_idx_type i = 0; i < M; ++i) {
                    B.xelem(i1[i] + 6 * (k * 3 + i2[i])) = H0d.xelem(k + N * i3[i]);
               }
          }
     }
private:
     const Material::MatType eMaterial;
     octave_idx_type iNumPreLoads;
     Matrix dTheta;
     NDArray epsilonRef;
     NDArray tauRef;
     ComplexMatrix f;
     Matrix e1, e2;
     const ElementData& oElemData;
};

class Iso8: public Element3D
{
public:
     Iso8(ElementTypes::TypeId eltype, octave_idx_type id, const Matrix& X, const Material* material, const int32NDArray& nodes, const ElementData& data)
          :Element3D(eltype, id, X, material, nodes, data) {
          FEM_ASSERT(nodes.numel() == 8);
     }

     static void AllocIntegrationRule(FemMatrixType eMatType) {
          static constexpr octave_idx_type N = 2;
          static constexpr double r[2][N] = {{0.577350269189626, -0.577350269189626}, {1., -1.}};
          static constexpr double alpha[2][N] = {{1., 1.}, {1., 1.}};

          const octave_idx_type iIntegRule = SelectIntegrationRule(eMatType);

          if (!rgIntegRule[iIntegRule].iGetNumEvalPoints()) {
               rgIntegRule[iIntegRule].SetNumEvalPoints(N * N * N, 3);

               octave_idx_type l = 0;

               for (octave_idx_type i = 0; i < N; ++i) {
                    for (octave_idx_type j = 0; j < N; ++j) {
                         for (octave_idx_type k = 0; k < N; ++k) {
                              rgIntegRule[iIntegRule].SetPosition(l, 0, r[iIntegRule][i]);
                              rgIntegRule[iIntegRule].SetPosition(l, 1, r[iIntegRule][j]);
                              rgIntegRule[iIntegRule].SetPosition(l, 2, r[iIntegRule][k]);
                              rgIntegRule[iIntegRule].SetWeight(l, alpha[iIntegRule][i] * alpha[iIntegRule][j] * alpha[iIntegRule][k]);
                              ++l;
                         }
                    }
               }
          }
     }

     virtual const IntegrationRule& GetIntegrationRule(FemMatrixType eMatType) const override final {
          const octave_idx_type iIntegRule = SelectIntegrationRule(eMatType);

          FEM_ASSERT(iIntegRule >= 0);
          FEM_ASSERT(static_cast<size_t>(iIntegRule) < rgIntegRule.size());
          FEM_ASSERT(rgIntegRule[iIntegRule].iGetNumEvalPoints() > 0);

          return rgIntegRule[iIntegRule];
     }

protected:
     virtual void ScalarInterpMatrixDer(const ColumnVector& rv, Matrix& Hd) const override final {
          const double r = rv.xelem(0);
          const double s = rv.xelem(1);
          const double t = rv.xelem(2);

          FEM_ASSERT(rv.numel() == 3);
          FEM_ASSERT(Hd.rows() == 8);
          FEM_ASSERT(Hd.columns() == 3);

          Hd.xelem(0) = ((s+1)*(t+1))/8.0E+0;
          Hd.xelem(1) = -((s+1)*(t+1))/8.0E+0;
          Hd.xelem(2) = -((1-s)*(t+1))/8.0E+0;
          Hd.xelem(3) = ((1-s)*(t+1))/8.0E+0;
          Hd.xelem(4) = ((s+1)*(1-t))/8.0E+0;
          Hd.xelem(5) = -((s+1)*(1-t))/8.0E+0;
          Hd.xelem(6) = -((1-s)*(1-t))/8.0E+0;
          Hd.xelem(7) = ((1-s)*(1-t))/8.0E+0;
          Hd.xelem(8) = ((r+1)*(t+1))/8.0E+0;
          Hd.xelem(9) = ((1-r)*(t+1))/8.0E+0;
          Hd.xelem(10) = -((1-r)*(t+1))/8.0E+0;
          Hd.xelem(11) = -((r+1)*(t+1))/8.0E+0;
          Hd.xelem(12) = ((r+1)*(1-t))/8.0E+0;
          Hd.xelem(13) = ((1-r)*(1-t))/8.0E+0;
          Hd.xelem(14) = -((1-r)*(1-t))/8.0E+0;
          Hd.xelem(15) = -((r+1)*(1-t))/8.0E+0;
          Hd.xelem(16) = ((r+1)*(s+1))/8.0E+0;
          Hd.xelem(17) = ((1-r)*(s+1))/8.0E+0;
          Hd.xelem(18) = ((1-r)*(1-s))/8.0E+0;
          Hd.xelem(19) = ((r+1)*(1-s))/8.0E+0;
          Hd.xelem(20) = -((r+1)*(s+1))/8.0E+0;
          Hd.xelem(21) = -((1-r)*(s+1))/8.0E+0;
          Hd.xelem(22) = -((1-r)*(1-s))/8.0E+0;
          Hd.xelem(23) = -((r+1)*(1-s))/8.0E+0;
     }

     virtual Matrix InterpGaussToNodal(FemMatrixType eMatType, const Matrix& taug) const override final {
          return InterpGaussToNodalTpl<double>(eMatType, taug);
     }

     virtual ComplexMatrix InterpGaussToNodal(FemMatrixType eMatType, const ComplexMatrix& taug) const override final {
          return InterpGaussToNodalTpl<std::complex<double> >(eMatType, taug);
     }

private:
     template <typename T>
     typename PostProcTypeTraits<T>::MatrixType
     InterpGaussToNodalTpl(FemMatrixType eMatType,
                           const typename PostProcTypeTraits<T>::MatrixType& taug) const {
          const IntegrationRule& oIntegRule = GetIntegrationRule(eMatType);
          const octave_idx_type iNumGauss = oIntegRule.iGetNumEvalPoints();
          const octave_idx_type iNumDir = oIntegRule.iGetNumDirections();
          const octave_idx_type iNumNodes = nodes.numel();

          ColumnVector rv(iNumDir);

          Matrix H(iNumGauss, iNumNodes);

          for (octave_idx_type i = 0; i < iNumGauss; ++i) {
               for (octave_idx_type j = 0; j < iNumDir; ++j) {
                    rv.xelem(j) = oIntegRule.dGetPosition(i, j);
               }

               ScalarInterpMatrix(rv, H, i);
          }

          return H.solve(taug);
     }

     void ScalarInterpMatrix(const ColumnVector& rv, Matrix& Hs, octave_idx_type irow) const override final {
          FEM_ASSERT(rv.numel() == 3);
          FEM_ASSERT(Hs.columns() == 8);
          FEM_ASSERT(irow >= 0);
          FEM_ASSERT(irow < Hs.rows());

          const double r = rv.xelem(0);
          const double s = rv.xelem(1);
          const double t = rv.xelem(2);
          const octave_idx_type nrows = Hs.rows();

          Hs.xelem(irow + nrows * 0) = ((r+1)*(s+1)*(t+1))/8.0;
          Hs.xelem(irow + nrows * 1) = ((1-r)*(s+1)*(t+1))/8.0;
          Hs.xelem(irow + nrows * 2) = ((1-r)*(1-s)*(t+1))/8.0;
          Hs.xelem(irow + nrows * 3) = ((r+1)*(1-s)*(t+1))/8.0;
          Hs.xelem(irow + nrows * 4) = ((r+1)*(s+1)*(1-t))/8.0;
          Hs.xelem(irow + nrows * 5) = ((1-r)*(s+1)*(1-t))/8.0;
          Hs.xelem(irow + nrows * 6) = ((1-r)*(1-s)*(1-t))/8.0;
          Hs.xelem(irow + nrows * 7) = ((r+1)*(1-s)*(1-t))/8.0;
     }

     static octave_idx_type SelectIntegrationRule(FemMatrixType eMatType) {
          switch (eMatType) {
          case MAT_MASS_LUMPED:
               return 1;
          default:
               return 0;
          }
     }

     static array<IntegrationRule, 2> rgIntegRule;
};

array<IntegrationRule, 2> Iso8::rgIntegRule;

class Iso20: public Element3D
{
public:
     Iso20(ElementTypes::TypeId eltype, octave_idx_type id, const Matrix& X, const Material* material, const int32NDArray& nodes, const ElementData& data)
          :Element3D(eltype, id, X, material, nodes, data) {
          FEM_ASSERT(nodes.numel() == 20);
     }

     static void AllocIntegrationRule(FemMatrixType eMatType) {
          constexpr octave_idx_type N = 3;
          static constexpr double r[2][N] = {{0.774596669241483, 0., -0.774596669241483}, {1., 0., -1.}};
          static constexpr double alpha[2][N] = {{0.555555555555556, 0.888888888888889, 0.555555555555556}, {2./3., 2./3., 2./3.}};

          const octave_idx_type iIntegRule = SelectIntegrationRule(eMatType);

          if (!rgIntegRule[iIntegRule].iGetNumEvalPoints()) {
               rgIntegRule[iIntegRule].SetNumEvalPoints(N * N * N, 3);

               octave_idx_type l = 0;

               for (octave_idx_type i = 0; i < N; ++i) {
                    for (octave_idx_type j = 0; j < N; ++j) {
                         for (octave_idx_type k = 0; k < N; ++k) {
                              rgIntegRule[iIntegRule].SetPosition(l, 0, r[iIntegRule][i]);
                              rgIntegRule[iIntegRule].SetPosition(l, 1, r[iIntegRule][j]);
                              rgIntegRule[iIntegRule].SetPosition(l, 2, r[iIntegRule][k]);
                              rgIntegRule[iIntegRule].SetWeight(l, alpha[iIntegRule][i] * alpha[iIntegRule][j] * alpha[iIntegRule][k]);
                              ++l;
                         }
                    }
               }
          }
     }

     virtual const IntegrationRule& GetIntegrationRule(FemMatrixType eMatType) const override final {
          const octave_idx_type iIntegRule = SelectIntegrationRule(eMatType);

          FEM_ASSERT(iIntegRule >= 0);
          FEM_ASSERT(static_cast<size_t>(iIntegRule) < rgIntegRule.size());
          FEM_ASSERT(rgIntegRule[iIntegRule].iGetNumEvalPoints() > 0);

          return rgIntegRule[iIntegRule];
     }

protected:
     virtual void ScalarInterpMatrixDer(const ColumnVector& rv, Matrix& Hd) const override final {
          const double r = rv.xelem(0);
          const double s = rv.xelem(1);
          const double t = rv.xelem(2);
          const double r2 = r * r;
          const double s2 = s * s;
          const double t2 = t * t;

          FEM_ASSERT(rv.numel() == 3);
          FEM_ASSERT(Hd.rows() == 20);
          FEM_ASSERT(Hd.columns() == 3);

          Hd.xelem(0) = ((s+1)*(t+1))/8.0E+0-(((s+1)*(1-t2))/4.0E+0+((1-s2)*(t+1))/4.0E+0-(r*(s+1)*(t+1))/2.0E+0)/2.0E+0;
          Hd.xelem(1) = (-((-((s+1)*(1-t2))/4.0E+0)-((1-s2)*(t+1))/4.0E+0-(r*(s+1)*(t+1))/2.0E+0)/2.0E+0)-((s+1)*(t+1))/8.0E+0;
          Hd.xelem(2) = (-((-((1-s)*(1-t2))/4.0E+0)-((1-s2)*(t+1))/4.0E+0-(r*(1-s)*(t+1))/2.0E+0)/2.0E+0)-((1-s)*(t+1))/8.0E+0;
          Hd.xelem(3) = ((1-s)*(t+1))/8.0E+0-(((1-s)*(1-t2))/4.0E+0+((1-s2)*(t+1))/4.0E+0-(r*(1-s)*(t+1))/2.0E+0)/2.0E+0;
          Hd.xelem(4) = ((s+1)*(1-t))/8.0E+0-(((s+1)*(1-t2))/4.0E+0+((1-s2)*(1-t))/4.0E+0-(r*(s+1)*(1-t))/2.0E+0)/2.0E+0;
          Hd.xelem(5) = (-((-((s+1)*(1-t2))/4.0E+0)-((1-s2)*(1-t))/4.0E+0-(r*(s+1)*(1-t))/2.0E+0)/2.0E+0)-((s+1)*(1-t))/8.0E+0;
          Hd.xelem(6) = (-((-((1-s)*(1-t2))/4.0E+0)-((1-s2)*(1-t))/4.0E+0-(r*(1-s)*(1-t))/2.0E+0)/2.0E+0)-((1-s)*(1-t))/8.0E+0;
          Hd.xelem(7) = ((1-s)*(1-t))/8.0E+0-(((1-s)*(1-t2))/4.0E+0+((1-s2)*(1-t))/4.0E+0-(r*(1-s)*(1-t))/2.0E+0)/2.0E+0;
          Hd.xelem(8) = -(r*(s+1)*(t+1))/2.0E+0;
          Hd.xelem(9) = -((1-s2)*(t+1))/4.0E+0;
          Hd.xelem(10) = -(r*(1-s)*(t+1))/2.0E+0;
          Hd.xelem(11) = ((1-s2)*(t+1))/4.0E+0;
          Hd.xelem(12) = -(r*(s+1)*(1-t))/2.0E+0;
          Hd.xelem(13) = -((1-s2)*(1-t))/4.0E+0;
          Hd.xelem(14) = -(r*(1-s)*(1-t))/2.0E+0;
          Hd.xelem(15) = ((1-s2)*(1-t))/4.0E+0;
          Hd.xelem(16) = ((s+1)*(1-t2))/4.0E+0;
          Hd.xelem(17) = -((s+1)*(1-t2))/4.0E+0;
          Hd.xelem(18) = -((1-s)*(1-t2))/4.0E+0;
          Hd.xelem(19) = ((1-s)*(1-t2))/4.0E+0;
          Hd.xelem(20) = ((r+1)*(t+1))/8.0E+0-(((r+1)*(1-t2))/4.0E+0-((r+1)*s*(t+1))/2.0E+0+((1-r2)*(t+1))/4.0E+0)/2.0E+0;
          Hd.xelem(21) = ((1-r)*(t+1))/8.0E+0-(((1-r)*(1-t2))/4.0E+0-((1-r)*s*(t+1))/2.0E+0+((1-r2)*(t+1))/4.0E+0)/2.0E+0;
          Hd.xelem(22) = (-((-((1-r)*(1-t2))/4.0E+0)-((1-r)*s*(t+1))/2.0E+0-((1-r2)*(t+1))/4.0E+0)/2.0E+0)-((1-r)*(t+1))/8.0E+0;
          Hd.xelem(23) = (-((-((r+1)*(1-t2))/4.0E+0)-((r+1)*s*(t+1))/2.0E+0-((1-r2)*(t+1))/4.0E+0)/2.0E+0)-((r+1)*(t+1))/8.0E+0;
          Hd.xelem(24) = ((r+1)*(1-t))/8.0E+0-(((r+1)*(1-t2))/4.0E+0-((r+1)*s*(1-t))/2.0E+0+((1-r2)*(1-t))/4.0E+0)/2.0E+0;
          Hd.xelem(25) = ((1-r)*(1-t))/8.0E+0-(((1-r)*(1-t2))/4.0E+0-((1-r)*s*(1-t))/2.0E+0+((1-r2)*(1-t))/4.0E+0)/2.0E+0;
          Hd.xelem(26) = (-((-((1-r)*(1-t2))/4.0E+0)-((1-r)*s*(1-t))/2.0E+0-((1-r2)*(1-t))/4.0E+0)/2.0E+0)-((1-r)*(1-t))/8.0E+0;
          Hd.xelem(27) = (-((-((r+1)*(1-t2))/4.0E+0)-((r+1)*s*(1-t))/2.0E+0-((1-r2)*(1-t))/4.0E+0)/2.0E+0)-((r+1)*(1-t))/8.0E+0;
          Hd.xelem(28) = ((1-r2)*(t+1))/4.0E+0;
          Hd.xelem(29) = -((1-r)*s*(t+1))/2.0E+0;
          Hd.xelem(30) = -((1-r2)*(t+1))/4.0E+0;
          Hd.xelem(31) = -((r+1)*s*(t+1))/2.0E+0;
          Hd.xelem(32) = ((1-r2)*(1-t))/4.0E+0;
          Hd.xelem(33) = -((1-r)*s*(1-t))/2.0E+0;
          Hd.xelem(34) = -((1-r2)*(1-t))/4.0E+0;
          Hd.xelem(35) = -((r+1)*s*(1-t))/2.0E+0;
          Hd.xelem(36) = ((r+1)*(1-t2))/4.0E+0;
          Hd.xelem(37) = ((1-r)*(1-t2))/4.0E+0;
          Hd.xelem(38) = -((1-r)*(1-t2))/4.0E+0;
          Hd.xelem(39) = -((r+1)*(1-t2))/4.0E+0;
          Hd.xelem(40) = ((r+1)*(s+1))/8.0E+0-((-((r+1)*(s+1)*t)/2.0E+0)+((r+1)*(1-s2))/4.0E+0+((1-r2)*(s+1))/4.0E+0)/2.0E+0;
          Hd.xelem(41) = ((1-r)*(s+1))/8.0E+0-((-((1-r)*(s+1)*t)/2.0E+0)+((1-r)*(1-s2))/4.0E+0+((1-r2)*(s+1))/4.0E+0)/2.0E+0;
          Hd.xelem(42) = ((1-r)*(1-s))/8.0E+0-((-((1-r)*(1-s)*t)/2.0E+0)+((1-r)*(1-s2))/4.0E+0+((1-r2)*(1-s))/4.0E+0)/2.0E+0;
          Hd.xelem(43) = ((r+1)*(1-s))/8.0E+0-((-((r+1)*(1-s)*t)/2.0E+0)+((r+1)*(1-s2))/4.0E+0+((1-r2)*(1-s))/4.0E+0)/2.0E+0;
          Hd.xelem(44) = (-((-((r+1)*(s+1)*t)/2.0E+0)-((r+1)*(1-s2))/4.0E+0-((1-r2)*(s+1))/4.0E+0)/2.0E+0)-((r+1)*(s+1))/8.0E+0;
          Hd.xelem(45) = (-((-((1-r)*(s+1)*t)/2.0E+0)-((1-r)*(1-s2))/4.0E+0-((1-r2)*(s+1))/4.0E+0)/2.0E+0)-((1-r)*(s+1))/8.0E+0;
          Hd.xelem(46) = (-((-((1-r)*(1-s)*t)/2.0E+0)-((1-r)*(1-s2))/4.0E+0-((1-r2)*(1-s))/4.0E+0)/2.0E+0)-((1-r)*(1-s))/8.0E+0;
          Hd.xelem(47) = (-((-((r+1)*(1-s)*t)/2.0E+0)-((r+1)*(1-s2))/4.0E+0-((1-r2)*(1-s))/4.0E+0)/2.0E+0)-((r+1)*(1-s))/8.0E+0;
          Hd.xelem(48) = ((1-r2)*(s+1))/4.0E+0;
          Hd.xelem(49) = ((1-r)*(1-s2))/4.0E+0;
          Hd.xelem(50) = ((1-r2)*(1-s))/4.0E+0;
          Hd.xelem(51) = ((r+1)*(1-s2))/4.0E+0;
          Hd.xelem(52) = -((1-r2)*(s+1))/4.0E+0;
          Hd.xelem(53) = -((1-r)*(1-s2))/4.0E+0;
          Hd.xelem(54) = -((1-r2)*(1-s))/4.0E+0;
          Hd.xelem(55) = -((r+1)*(1-s2))/4.0E+0;
          Hd.xelem(56) = -((r+1)*(s+1)*t)/2.0E+0;
          Hd.xelem(57) = -((1-r)*(s+1)*t)/2.0E+0;
          Hd.xelem(58) = -((1-r)*(1-s)*t)/2.0E+0;
          Hd.xelem(59) = -((r+1)*(1-s)*t)/2.0E+0;
     }

     virtual Matrix InterpGaussToNodal(FemMatrixType eMatType, const Matrix& taug) const override final {
          return InterpGaussToNodalTpl<double>(eMatType, taug);
     }

     virtual ComplexMatrix InterpGaussToNodal(FemMatrixType eMatType, const ComplexMatrix& taug) const override final {
          return InterpGaussToNodalTpl<std::complex<double> >(eMatType, taug);
     }

private:
     template <typename T>
     typename PostProcTypeTraits<T>::MatrixType
     InterpGaussToNodalTpl(FemMatrixType eMatType,
                           const typename PostProcTypeTraits<T>::MatrixType& taug) const {
          const IntegrationRule& oIntegRule = GetIntegrationRule(eMatType);
          const octave_idx_type iNumGauss = oIntegRule.iGetNumEvalPoints();
          const octave_idx_type iNumDir = oIntegRule.iGetNumDirections();
          const octave_idx_type iNumNodes = nodes.numel();
          ColumnVector rv(iNumDir);

          Matrix H(iNumGauss, iNumNodes);

          for (octave_idx_type i = 0; i < iNumGauss; ++i) {
               for (octave_idx_type j = 0; j < iNumDir; ++j) {
                    rv.xelem(j) = oIntegRule.dGetPosition(i, j);
               }

               ScalarInterpMatrix(rv, H, i);
          }

          return H.solve(taug);
     }

     void ScalarInterpMatrix(const ColumnVector& rv, Matrix& Hs, octave_idx_type irow) const override final {
          FEM_ASSERT(rv.numel() == 3);
          FEM_ASSERT(Hs.columns() == 20);
          FEM_ASSERT(irow >= 0);
          FEM_ASSERT(irow < Hs.rows());

          const double r = rv.xelem(0);
          const double s = rv.xelem(1);
          const double t = rv.xelem(2);
          const double r2 = r * r;
          const double s2 = s * s;
          const double t2 = t * t;
          const octave_idx_type nrows = Hs.rows();

          Hs.xelem(irow + nrows * 0) = ((r+1)*(s+1)*(t+1))/8.0E+0-(((r+1)*(s+1)*(1-t2))/4.0E+0+((r+1)*(1-s2)*(t+1))/4.0E+0+((1-r2)*(s+1)*(t+1))/4.0E+0)/2.0E+0;
          Hs.xelem(irow + nrows * 1) = ((1-r)*(s+1)*(t+1))/8.0E+0-(((1-r)*(s+1)*(1-t2))/4.0E+0+((1-r)*(1-s2)*(t+1))/4.0E+0+((1-r2)*(s+1)*(t+1))/4.0E+0)/2.0E+0;
          Hs.xelem(irow + nrows * 2) = ((1-r)*(1-s)*(t+1))/8.0E+0-(((1-r)*(1-s)*(1-t2))/4.0E+0+((1-r)*(1-s2)*(t+1))/4.0E+0+((1-r2)*(1-s)*(t+1))/4.0E+0)/2.0E+0;
          Hs.xelem(irow + nrows * 3) = ((r+1)*(1-s)*(t+1))/8.0E+0-(((r+1)*(1-s)*(1-t2))/4.0E+0+((r+1)*(1-s2)*(t+1))/4.0E+0+((1-r2)*(1-s)*(t+1))/4.0E+0)/2.0E+0;
          Hs.xelem(irow + nrows * 4) = ((r+1)*(s+1)*(1-t))/8.0E+0-(((r+1)*(s+1)*(1-t2))/4.0E+0+((r+1)*(1-s2)*(1-t))/4.0E+0+((1-r2)*(s+1)*(1-t))/4.0E+0)/2.0E+0;
          Hs.xelem(irow + nrows * 5) = ((1-r)*(s+1)*(1-t))/8.0E+0-(((1-r)*(s+1)*(1-t2))/4.0E+0+((1-r)*(1-s2)*(1-t))/4.0E+0+((1-r2)*(s+1)*(1-t))/4.0E+0)/2.0E+0;
          Hs.xelem(irow + nrows * 6) = ((1-r)*(1-s)*(1-t))/8.0E+0-(((1-r)*(1-s)*(1-t2))/4.0E+0+((1-r)*(1-s2)*(1-t))/4.0E+0+((1-r2)*(1-s)*(1-t))/4.0E+0)/2.0E+0;
          Hs.xelem(irow + nrows * 7) = ((r+1)*(1-s)*(1-t))/8.0E+0-(((r+1)*(1-s)*(1-t2))/4.0E+0+((r+1)*(1-s2)*(1-t))/4.0E+0+((1-r2)*(1-s)*(1-t))/4.0E+0)/2.0E+0;
          Hs.xelem(irow + nrows * 8) = ((1-r2)*(s+1)*(t+1))/4.0E+0;
          Hs.xelem(irow + nrows * 9) = ((1-r)*(1-s2)*(t+1))/4.0E+0;
          Hs.xelem(irow + nrows * 10) = ((1-r2)*(1-s)*(t+1))/4.0E+0;
          Hs.xelem(irow + nrows * 11) = ((r+1)*(1-s2)*(t+1))/4.0E+0;
          Hs.xelem(irow + nrows * 12) = ((1-r2)*(s+1)*(1-t))/4.0E+0;
          Hs.xelem(irow + nrows * 13) = ((1-r)*(1-s2)*(1-t))/4.0E+0;
          Hs.xelem(irow + nrows * 14) = ((1-r2)*(1-s)*(1-t))/4.0E+0;
          Hs.xelem(irow + nrows * 15) = ((r+1)*(1-s2)*(1-t))/4.0E+0;
          Hs.xelem(irow + nrows * 16) = ((r+1)*(s+1)*(1-t2))/4.0E+0;
          Hs.xelem(irow + nrows * 17) = ((1-r)*(s+1)*(1-t2))/4.0E+0;
          Hs.xelem(irow + nrows * 18) = ((1-r)*(1-s)*(1-t2))/4.0E+0;
          Hs.xelem(irow + nrows * 19) = ((r+1)*(1-s)*(1-t2))/4.0E+0;
     }

     static octave_idx_type SelectIntegrationRule(FemMatrixType eMatType) {
          switch (eMatType) {
          case MAT_MASS_LUMPED:
               return 1;
          default:
               return 0;
          }
     }

     static array<IntegrationRule, 2> rgIntegRule;
};

array<IntegrationRule, 2> Iso20::rgIntegRule;

class Iso20r: public Element3D
{
public:
     Iso20r(ElementTypes::TypeId eltype, octave_idx_type id, const Matrix& X, const Material* material, const int32NDArray& nodes, const ElementData& data)
          :Element3D(eltype, id, X, material, nodes, data) {
          FEM_ASSERT(nodes.numel() == 20);
     }

     static void AllocIntegrationRule(FemMatrixType eMatType) {
          constexpr octave_idx_type NMAX = 3;
          constexpr octave_idx_type NINTEG = 3;
          static constexpr octave_idx_type N[NINTEG] = {2, 3, 3};
          static constexpr double r[NINTEG][NMAX] = {{0.577350269189626, -0.577350269189626}, {0.774596669241483, 0., -0.774596669241483}, {1., 0., -1.}};
          static constexpr double alpha[NINTEG][NMAX] = {{1., 1.}, {0.555555555555556, 0.888888888888889, 0.555555555555556}, {2./3., 2./3., 2./3.}};

          const octave_idx_type iIntegRule = SelectIntegrationRule(eMatType);

          if (!rgIntegRule[iIntegRule].iGetNumEvalPoints()) {
               rgIntegRule[iIntegRule].SetNumEvalPoints(N[iIntegRule] * N[iIntegRule] * N[iIntegRule], 3);

               octave_idx_type l = 0;

               for (octave_idx_type i = 0; i < N[iIntegRule]; ++i) {
                    for (octave_idx_type j = 0; j < N[iIntegRule]; ++j) {
                         for (octave_idx_type k = 0; k < N[iIntegRule]; ++k) {
                              rgIntegRule[iIntegRule].SetPosition(l, 0, r[iIntegRule][i]);
                              rgIntegRule[iIntegRule].SetPosition(l, 1, r[iIntegRule][j]);
                              rgIntegRule[iIntegRule].SetPosition(l, 2, r[iIntegRule][k]);
                              rgIntegRule[iIntegRule].SetWeight(l, alpha[iIntegRule][i] * alpha[iIntegRule][j] * alpha[iIntegRule][k]);
                              ++l;
                         }
                    }
               }
          }
     }

     virtual const IntegrationRule& GetIntegrationRule(FemMatrixType eMatType) const override final {
          const octave_idx_type iIntegRule = SelectIntegrationRule(eMatType);

          FEM_ASSERT(iIntegRule >= 0);
          FEM_ASSERT(static_cast<size_t>(iIntegRule) < rgIntegRule.size());
          FEM_ASSERT(rgIntegRule[iIntegRule].iGetNumEvalPoints() > 0);

          return rgIntegRule[iIntegRule];
     }

protected:
     virtual void ScalarInterpMatrixDer(const ColumnVector& rv, Matrix& Hd) const override final {
          const double r1 = rv.xelem(0);
          const double r2 = rv.xelem(1);
          const double r3 = rv.xelem(2);

          FEM_ASSERT(rv.numel() == 3);
          FEM_ASSERT(Hd.rows() == 20);
          FEM_ASSERT(Hd.columns() == 3);

          Hd.xelem(0) = 1.25E-1*(1.0E+0-r2)*(1.0E+0-r3)*(r3+r2+r1+2.0E+0)+1.25E-1*(r1-1.0E+0)*(1.0E+0-r2)*(1.0E+0-r3);
          Hd.xelem(1) = (-1.25E-1*(1.0E+0-r2)*(1.0E+0-r3)*(r3+r2-r1+2.0E+0))-1.25E-1*((-r1)-1.0E+0)*(1.0E+0-r2)*(1.0E+0-r3);
          Hd.xelem(2) = (-1.25E-1*(r2+1.0E+0)*(1.0E+0-r3)*(r3-r2-r1+2.0E+0))-1.25E-1*((-r1)-1.0E+0)*(r2+1.0E+0)*(1.0E+0-r3);
          Hd.xelem(3) = 1.25E-1*(r2+1.0E+0)*(1.0E+0-r3)*(r3-r2+r1+2.0E+0)+1.25E-1*(r1-1.0E+0)*(r2+1.0E+0)*(1.0E+0-r3);
          Hd.xelem(4) = 1.25E-1*(1.0E+0-r2)*((-r3)+r2+r1+2.0E+0)*(r3+1.0E+0)+1.25E-1*(r1-1.0E+0)*(1.0E+0-r2)*(r3+1.0E+0);
          Hd.xelem(5) = (-1.25E-1*(1.0E+0-r2)*((-r3)+r2-r1+2.0E+0)*(r3+1.0E+0))-1.25E-1*((-r1)-1.0E+0)*(1.0E+0-r2)*(r3+1.0E+0);
          Hd.xelem(6) = (-1.25E-1*(r2+1.0E+0)*((-r3)-r2-r1+2.0E+0)*(r3+1.0E+0))-1.25E-1*((-r1)-1.0E+0)*(r2+1.0E+0)*(r3+1.0E+0);
          Hd.xelem(7) = 1.25E-1*(r2+1.0E+0)*((-r3)-r2+r1+2.0E+0)*(r3+1.0E+0)+1.25E-1*(r1-1.0E+0)*(r2+1.0E+0)*(r3+1.0E+0);
          Hd.xelem(8) = 2.5E-1*(1.0E+0-r1)*(1.0E+0-r2)*(1.0E+0-r3)-2.5E-1*(r1+1.0E+0)*(1.0E+0-r2)*(1.0E+0-r3);
          Hd.xelem(9) = 2.5E-1*(1.0E+0-r2)*(r2+1.0E+0)*(1.0E+0-r3);
          Hd.xelem(10) = 2.5E-1*(1.0E+0-r1)*(r2+1.0E+0)*(1.0E+0-r3)-2.5E-1*(r1+1.0E+0)*(r2+1.0E+0)*(1.0E+0-r3);
          Hd.xelem(11) = -2.5E-1*(1.0E+0-r2)*(r2+1.0E+0)*(1.0E+0-r3);
          Hd.xelem(12) = 2.5E-1*(1.0E+0-r1)*(1.0E+0-r2)*(r3+1.0E+0)-2.5E-1*(r1+1.0E+0)*(1.0E+0-r2)*(r3+1.0E+0);
          Hd.xelem(13) = 2.5E-1*(1.0E+0-r2)*(r2+1.0E+0)*(r3+1.0E+0);
          Hd.xelem(14) = 2.5E-1*(1.0E+0-r1)*(r2+1.0E+0)*(r3+1.0E+0)-2.5E-1*(r1+1.0E+0)*(r2+1.0E+0)*(r3+1.0E+0);
          Hd.xelem(15) = -2.5E-1*(1.0E+0-r2)*(r2+1.0E+0)*(r3+1.0E+0);
          Hd.xelem(16) = -2.5E-1*(1.0E+0-r2)*(1.0E+0-r3)*(r3+1.0E+0);
          Hd.xelem(17) = 2.5E-1*(1.0E+0-r2)*(1.0E+0-r3)*(r3+1.0E+0);
          Hd.xelem(18) = 2.5E-1*(r2+1.0E+0)*(1.0E+0-r3)*(r3+1.0E+0);
          Hd.xelem(19) = -2.5E-1*(r2+1.0E+0)*(1.0E+0-r3)*(r3+1.0E+0);
          Hd.xelem(20) = 1.25E-1*(r1-1.0E+0)*(1.0E+0-r2)*(1.0E+0-r3)-1.25E-1*(r1-1.0E+0)*(1.0E+0-r3)*(r3+r2+r1+2.0E+0);
          Hd.xelem(21) = 1.25E-1*((-r1)-1.0E+0)*(1.0E+0-r2)*(1.0E+0-r3)-1.25E-1*((-r1)-1.0E+0)*(1.0E+0-r3)*(r3+r2-r1+2.0E+0);
          Hd.xelem(22) = 1.25E-1*((-r1)-1.0E+0)*(1.0E+0-r3)*(r3-r2-r1+2.0E+0)-1.25E-1*((-r1)-1.0E+0)*(r2+1.0E+0)*(1.0E+0-r3);
          Hd.xelem(23) = 1.25E-1*(r1-1.0E+0)*(1.0E+0-r3)*(r3-r2+r1+2.0E+0)-1.25E-1*(r1-1.0E+0)*(r2+1.0E+0)*(1.0E+0-r3);
          Hd.xelem(24) = 1.25E-1*(r1-1.0E+0)*(1.0E+0-r2)*(r3+1.0E+0)-1.25E-1*(r1-1.0E+0)*((-r3)+r2+r1+2.0E+0)*(r3+1.0E+0);
          Hd.xelem(25) = 1.25E-1*((-r1)-1.0E+0)*(1.0E+0-r2)*(r3+1.0E+0)-1.25E-1*((-r1)-1.0E+0)*((-r3)+r2-r1+2.0E+0)*(r3+1.0E+0);
          Hd.xelem(26) = 1.25E-1*((-r1)-1.0E+0)*((-r3)-r2-r1+2.0E+0)*(r3+1.0E+0)-1.25E-1*((-r1)-1.0E+0)*(r2+1.0E+0)*(r3+1.0E+0);
          Hd.xelem(27) = 1.25E-1*(r1-1.0E+0)*((-r3)-r2+r1+2.0E+0)*(r3+1.0E+0)-1.25E-1*(r1-1.0E+0)*(r2+1.0E+0)*(r3+1.0E+0);
          Hd.xelem(28) = -2.5E-1*(1.0E+0-r1)*(r1+1.0E+0)*(1.0E+0-r3);
          Hd.xelem(29) = 2.5E-1*(r1+1.0E+0)*(1.0E+0-r2)*(1.0E+0-r3)-2.5E-1*(r1+1.0E+0)*(r2+1.0E+0)*(1.0E+0-r3);
          Hd.xelem(30) = 2.5E-1*(1.0E+0-r1)*(r1+1.0E+0)*(1.0E+0-r3);
          Hd.xelem(31) = 2.5E-1*(1.0E+0-r1)*(1.0E+0-r2)*(1.0E+0-r3)-2.5E-1*(1.0E+0-r1)*(r2+1.0E+0)*(1.0E+0-r3);
          Hd.xelem(32) = -2.5E-1*(1.0E+0-r1)*(r1+1.0E+0)*(r3+1.0E+0);
          Hd.xelem(33) = 2.5E-1*(r1+1.0E+0)*(1.0E+0-r2)*(r3+1.0E+0)-2.5E-1*(r1+1.0E+0)*(r2+1.0E+0)*(r3+1.0E+0);
          Hd.xelem(34) = 2.5E-1*(1.0E+0-r1)*(r1+1.0E+0)*(r3+1.0E+0);
          Hd.xelem(35) = 2.5E-1*(1.0E+0-r1)*(1.0E+0-r2)*(r3+1.0E+0)-2.5E-1*(1.0E+0-r1)*(r2+1.0E+0)*(r3+1.0E+0);
          Hd.xelem(36) = -2.5E-1*(1.0E+0-r1)*(1.0E+0-r3)*(r3+1.0E+0);
          Hd.xelem(37) = -2.5E-1*(r1+1.0E+0)*(1.0E+0-r3)*(r3+1.0E+0);
          Hd.xelem(38) = 2.5E-1*(r1+1.0E+0)*(1.0E+0-r3)*(r3+1.0E+0);
          Hd.xelem(39) = 2.5E-1*(1.0E+0-r1)*(1.0E+0-r3)*(r3+1.0E+0);
          Hd.xelem(40) = 1.25E-1*(r1-1.0E+0)*(1.0E+0-r2)*(1.0E+0-r3)-1.25E-1*(r1-1.0E+0)*(1.0E+0-r2)*(r3+r2+r1+2.0E+0);
          Hd.xelem(41) = 1.25E-1*((-r1)-1.0E+0)*(1.0E+0-r2)*(1.0E+0-r3)-1.25E-1*((-r1)-1.0E+0)*(1.0E+0-r2)*(r3+r2-r1+2.0E+0);
          Hd.xelem(42) = 1.25E-1*((-r1)-1.0E+0)*(r2+1.0E+0)*(1.0E+0-r3)-1.25E-1*((-r1)-1.0E+0)*(r2+1.0E+0)*(r3-r2-r1+2.0E+0);
          Hd.xelem(43) = 1.25E-1*(r1-1.0E+0)*(r2+1.0E+0)*(1.0E+0-r3)-1.25E-1*(r1-1.0E+0)*(r2+1.0E+0)*(r3-r2+r1+2.0E+0);
          Hd.xelem(44) = 1.25E-1*(r1-1.0E+0)*(1.0E+0-r2)*((-r3)+r2+r1+2.0E+0)-1.25E-1*(r1-1.0E+0)*(1.0E+0-r2)*(r3+1.0E+0);
          Hd.xelem(45) = 1.25E-1*((-r1)-1.0E+0)*(1.0E+0-r2)*((-r3)+r2-r1+2.0E+0)-1.25E-1*((-r1)-1.0E+0)*(1.0E+0-r2)*(r3+1.0E+0);
          Hd.xelem(46) = 1.25E-1*((-r1)-1.0E+0)*(r2+1.0E+0)*((-r3)-r2-r1+2.0E+0)-1.25E-1*((-r1)-1.0E+0)*(r2+1.0E+0)*(r3+1.0E+0);
          Hd.xelem(47) = 1.25E-1*(r1-1.0E+0)*(r2+1.0E+0)*((-r3)-r2+r1+2.0E+0)-1.25E-1*(r1-1.0E+0)*(r2+1.0E+0)*(r3+1.0E+0);
          Hd.xelem(48) = -2.5E-1*(1.0E+0-r1)*(r1+1.0E+0)*(1.0E+0-r2);
          Hd.xelem(49) = -2.5E-1*(r1+1.0E+0)*(1.0E+0-r2)*(r2+1.0E+0);
          Hd.xelem(50) = -2.5E-1*(1.0E+0-r1)*(r1+1.0E+0)*(r2+1.0E+0);
          Hd.xelem(51) = -2.5E-1*(1.0E+0-r1)*(1.0E+0-r2)*(r2+1.0E+0);
          Hd.xelem(52) = 2.5E-1*(1.0E+0-r1)*(r1+1.0E+0)*(1.0E+0-r2);
          Hd.xelem(53) = 2.5E-1*(r1+1.0E+0)*(1.0E+0-r2)*(r2+1.0E+0);
          Hd.xelem(54) = 2.5E-1*(1.0E+0-r1)*(r1+1.0E+0)*(r2+1.0E+0);
          Hd.xelem(55) = 2.5E-1*(1.0E+0-r1)*(1.0E+0-r2)*(r2+1.0E+0);
          Hd.xelem(56) = 2.5E-1*(1.0E+0-r1)*(1.0E+0-r2)*(1.0E+0-r3)-2.5E-1*(1.0E+0-r1)*(1.0E+0-r2)*(r3+1.0E+0);
          Hd.xelem(57) = 2.5E-1*(r1+1.0E+0)*(1.0E+0-r2)*(1.0E+0-r3)-2.5E-1*(r1+1.0E+0)*(1.0E+0-r2)*(r3+1.0E+0);
          Hd.xelem(58) = 2.5E-1*(r1+1.0E+0)*(r2+1.0E+0)*(1.0E+0-r3)-2.5E-1*(r1+1.0E+0)*(r2+1.0E+0)*(r3+1.0E+0);
          Hd.xelem(59) = 2.5E-1*(1.0E+0-r1)*(r2+1.0E+0)*(1.0E+0-r3)-2.5E-1*(1.0E+0-r1)*(r2+1.0E+0)*(r3+1.0E+0);
     }

     virtual Matrix InterpGaussToNodal(FemMatrixType eMatType, const Matrix& taug) const override final {
          return InterpGaussToNodalTpl<double>(eMatType, taug);
     }

     virtual ComplexMatrix InterpGaussToNodal(FemMatrixType eMatType, const ComplexMatrix& taug) const override final {
          return InterpGaussToNodalTpl<std::complex<double> >(eMatType, taug);
     }

private:
     template <typename T>
     typename PostProcTypeTraits<T>::MatrixType
     InterpGaussToNodalTpl(FemMatrixType eMatType,
                           const typename PostProcTypeTraits<T>::MatrixType& taug) const {
          const IntegrationRule& oIntegRule = GetIntegrationRule(eMatType);
          const octave_idx_type iNumGauss = oIntegRule.iGetNumEvalPoints();
          const octave_idx_type iNumDir = oIntegRule.iGetNumDirections();
          const octave_idx_type iNumNodes = nodes.numel();
          ColumnVector rv(iNumDir);

          Matrix H(iNumGauss, iNumNodes);

          for (octave_idx_type i = 0; i < iNumGauss; ++i) {
               for (octave_idx_type j = 0; j < iNumDir; ++j) {
                    rv.xelem(j) = oIntegRule.dGetPosition(i, j);
               }

               ScalarInterpMatrix(rv, H, i);
          }

          return H.solve(taug);
     }

     void ScalarInterpMatrix(const ColumnVector& rv, Matrix& Hs, octave_idx_type irow) const override final {
          FEM_ASSERT(rv.numel() == 3);
          FEM_ASSERT(Hs.columns() == 20);
          FEM_ASSERT(irow >= 0);
          FEM_ASSERT(irow < Hs.rows());

          const double r1 = rv.xelem(0);
          const double r2 = rv.xelem(1);
          const double r3 = rv.xelem(2);
          const octave_idx_type nrows = Hs.rows();

          Hs.xelem(irow + nrows * 0) = 1.25E-1*(r1-1.0E+0)*(1.0E+0-r2)*(1.0E+0-r3)*(r3+r2+r1+2.0E+0);
          Hs.xelem(irow + nrows * 1) = 1.25E-1*((-r1)-1.0E+0)*(1.0E+0-r2)*(1.0E+0-r3)*(r3+r2-r1+2.0E+0);
          Hs.xelem(irow + nrows * 2) = 1.25E-1*((-r1)-1.0E+0)*(r2+1.0E+0)*(1.0E+0-r3)*(r3-r2-r1+2.0E+0);
          Hs.xelem(irow + nrows * 3) = 1.25E-1*(r1-1.0E+0)*(r2+1.0E+0)*(1.0E+0-r3)*(r3-r2+r1+2.0E+0);
          Hs.xelem(irow + nrows * 4) = 1.25E-1*(r1-1.0E+0)*(1.0E+0-r2)*((-r3)+r2+r1+2.0E+0)*(r3+1.0E+0);
          Hs.xelem(irow + nrows * 5) = 1.25E-1*((-r1)-1.0E+0)*(1.0E+0-r2)*((-r3)+r2-r1+2.0E+0)*(r3+1.0E+0);
          Hs.xelem(irow + nrows * 6) = 1.25E-1*((-r1)-1.0E+0)*(r2+1.0E+0)*((-r3)-r2-r1+2.0E+0)*(r3+1.0E+0);
          Hs.xelem(irow + nrows * 7) = 1.25E-1*(r1-1.0E+0)*(r2+1.0E+0)*((-r3)-r2+r1+2.0E+0)*(r3+1.0E+0);
          Hs.xelem(irow + nrows * 8) = 2.5E-1*(1.0E+0-r1)*(r1+1.0E+0)*(1.0E+0-r2)*(1.0E+0-r3);
          Hs.xelem(irow + nrows * 9) = 2.5E-1*(r1+1.0E+0)*(1.0E+0-r2)*(r2+1.0E+0)*(1.0E+0-r3);
          Hs.xelem(irow + nrows * 10) = 2.5E-1*(1.0E+0-r1)*(r1+1.0E+0)*(r2+1.0E+0)*(1.0E+0-r3);
          Hs.xelem(irow + nrows * 11) = 2.5E-1*(1.0E+0-r1)*(1.0E+0-r2)*(r2+1.0E+0)*(1.0E+0-r3);
          Hs.xelem(irow + nrows * 12) = 2.5E-1*(1.0E+0-r1)*(r1+1.0E+0)*(1.0E+0-r2)*(r3+1.0E+0);
          Hs.xelem(irow + nrows * 13) = 2.5E-1*(r1+1.0E+0)*(1.0E+0-r2)*(r2+1.0E+0)*(r3+1.0E+0);
          Hs.xelem(irow + nrows * 14) = 2.5E-1*(1.0E+0-r1)*(r1+1.0E+0)*(r2+1.0E+0)*(r3+1.0E+0);
          Hs.xelem(irow + nrows * 15) = 2.5E-1*(1.0E+0-r1)*(1.0E+0-r2)*(r2+1.0E+0)*(r3+1.0E+0);
          Hs.xelem(irow + nrows * 16) = 2.5E-1*(1.0E+0-r1)*(1.0E+0-r2)*(1.0E+0-r3)*(r3+1.0E+0);
          Hs.xelem(irow + nrows * 17) = 2.5E-1*(r1+1.0E+0)*(1.0E+0-r2)*(1.0E+0-r3)*(r3+1.0E+0);
          Hs.xelem(irow + nrows * 18) = 2.5E-1*(r1+1.0E+0)*(r2+1.0E+0)*(1.0E+0-r3)*(r3+1.0E+0);
          Hs.xelem(irow + nrows * 19) = 2.5E-1*(1.0E+0-r1)*(r2+1.0E+0)*(1.0E+0-r3)*(r3+1.0E+0);

     }

     static octave_idx_type SelectIntegrationRule(FemMatrixType eMatType) {
          switch (eMatType) {
          case MAT_MASS_LUMPED:
               return 2;
          case MAT_STIFFNESS:
          case MAT_STIFFNESS_IM:
          case MAT_STIFFNESS_SYM:
          case MAT_STIFFNESS_SYM_L:
          case MAT_STIFFNESS_TAU0:
          case MAT_STIFFNESS_OMEGA:
          case MAT_STIFFNESS_OMEGA_DOT:
          case MAT_STIFFNESS_ACOUSTICS_RE:
          case MAT_STIFFNESS_ACOUSTICS_IM:
          case MAT_STIFFNESS_FLUID_STRUCT_RE:
          case MAT_STIFFNESS_FLUID_STRUCT_IM:
          case VEC_COLL_STIFFNESS:
          case VEC_COLL_STIFF_ACOUSTICS:
          case MAT_DAMPING:
          case MAT_DAMPING_SYM:
          case MAT_DAMPING_SYM_L:
          case MAT_DAMPING_ACOUSTICS_RE:
          case MAT_DAMPING_ACOUSTICS_IM:
          case MAT_DAMPING_FLUID_STRUCT_RE:
          case MAT_DAMPING_FLUID_STRUCT_IM:
          case MAT_DAMPING_OMEGA:
               return 0;
          default:
               return 1;
          }
     }

     static array<IntegrationRule, 3> rgIntegRule;
};

array<IntegrationRule, 3> Iso20r::rgIntegRule;

class Iso27: public Element3D
{
public:
     Iso27(ElementTypes::TypeId eltype, octave_idx_type id, const Matrix& X, const Material* material, const int32NDArray& nodes, const ElementData& data)
          :Element3D(eltype, id, X, material, nodes, data) {
          FEM_ASSERT(nodes.numel() == 27);
     }

     static void AllocIntegrationRule(FemMatrixType eMatType) {
          constexpr octave_idx_type N = 3;
          static constexpr double r[2][N] = {{0.774596669241483, 0., -0.774596669241483}, {1., 0., -1.}};
          static constexpr double alpha[2][N] = {{0.555555555555556, 0.888888888888889, 0.555555555555556}, {2./3., 2./3., 2./3.}};

          const octave_idx_type iIntegRule = SelectIntegrationRule(eMatType);

          if (!rgIntegRule[iIntegRule].iGetNumEvalPoints()) {
               rgIntegRule[iIntegRule].SetNumEvalPoints(N * N * N, 3);

               octave_idx_type l = 0;

               for (octave_idx_type i = 0; i < N; ++i) {
                    for (octave_idx_type j = 0; j < N; ++j) {
                         for (octave_idx_type k = 0; k < N; ++k) {
                              rgIntegRule[iIntegRule].SetPosition(l, 0, r[iIntegRule][i]);
                              rgIntegRule[iIntegRule].SetPosition(l, 1, r[iIntegRule][j]);
                              rgIntegRule[iIntegRule].SetPosition(l, 2, r[iIntegRule][k]);
                              rgIntegRule[iIntegRule].SetWeight(l, alpha[iIntegRule][i] * alpha[iIntegRule][j] * alpha[iIntegRule][k]);
                              ++l;
                         }
                    }
               }
          }
     }

     virtual const IntegrationRule& GetIntegrationRule(FemMatrixType eMatType) const override final {
          const octave_idx_type iIntegRule = SelectIntegrationRule(eMatType);

          FEM_ASSERT(iIntegRule >= 0);
          FEM_ASSERT(static_cast<size_t>(iIntegRule) < rgIntegRule.size());
          FEM_ASSERT(rgIntegRule[iIntegRule].iGetNumEvalPoints() > 0);

          return rgIntegRule[iIntegRule];
     }

protected:
     virtual void ScalarInterpMatrixDer(const ColumnVector& rv, Matrix& Hd) const override final {
          const double r = rv.xelem(0);
          const double s = rv.xelem(1);
          const double t = rv.xelem(2);
          const double r2 = r * r;
          const double s2 = s * s;
          const double t2 = t * t;

          FEM_ASSERT(rv.numel() == 3);
          FEM_ASSERT(Hd.rows() == 27);
          FEM_ASSERT(Hd.columns() == 3);

          Hd.xelem(0) = (r*(s-1)*s*(t-1)*t)/8.0E+0+((r-1)*(s-1)*s*(t-1)*t)/8.0E+0;
          Hd.xelem(1) = ((r+1)*(s-1)*s*(t-1)*t)/8.0E+0+(r*(s-1)*s*(t-1)*t)/8.0E+0;
          Hd.xelem(2) = ((r+1)*s*(s+1)*(t-1)*t)/8.0E+0+(r*s*(s+1)*(t-1)*t)/8.0E+0;
          Hd.xelem(3) = (r*s*(s+1)*(t-1)*t)/8.0E+0+((r-1)*s*(s+1)*(t-1)*t)/8.0E+0;
          Hd.xelem(4) = (r*(s-1)*s*t*(t+1))/8.0E+0+((r-1)*(s-1)*s*t*(t+1))/8.0E+0;
          Hd.xelem(5) = ((r+1)*(s-1)*s*t*(t+1))/8.0E+0+(r*(s-1)*s*t*(t+1))/8.0E+0;
          Hd.xelem(6) = ((r+1)*s*(s+1)*t*(t+1))/8.0E+0+(r*s*(s+1)*t*(t+1))/8.0E+0;
          Hd.xelem(7) = (r*s*(s+1)*t*(t+1))/8.0E+0+((r-1)*s*(s+1)*t*(t+1))/8.0E+0;
          Hd.xelem(8) = -(r*(s-1)*s*(t-1)*t)/2.0E+0;
          Hd.xelem(9) = ((r+1)*(1-s2)*(t-1)*t)/4.0E+0+(r*(1-s2)*(t-1)*t)/4.0E+0;
          Hd.xelem(10) = -(r*s*(s+1)*(t-1)*t)/2.0E+0;
          Hd.xelem(11) = (r*(1-s2)*(t-1)*t)/4.0E+0+((r-1)*(1-s2)*(t-1)*t)/4.0E+0;
          Hd.xelem(12) = (r*(s-1)*s*(1-t2))/4.0E+0+((r-1)*(s-1)*s*(1-t2))/4.0E+0;
          Hd.xelem(13) = ((r+1)*(s-1)*s*(1-t2))/4.0E+0+(r*(s-1)*s*(1-t2))/4.0E+0;
          Hd.xelem(14) = ((r+1)*s*(s+1)*(1-t2))/4.0E+0+(r*s*(s+1)*(1-t2))/4.0E+0;
          Hd.xelem(15) = (r*s*(s+1)*(1-t2))/4.0E+0+((r-1)*s*(s+1)*(1-t2))/4.0E+0;
          Hd.xelem(16) = -(r*(s-1)*s*t*(t+1))/2.0E+0;
          Hd.xelem(17) = ((r+1)*(1-s2)*t*(t+1))/4.0E+0+(r*(1-s2)*t*(t+1))/4.0E+0;
          Hd.xelem(18) = -(r*s*(s+1)*t*(t+1))/2.0E+0;
          Hd.xelem(19) = (r*(1-s2)*t*(t+1))/4.0E+0+((r-1)*(1-s2)*t*(t+1))/4.0E+0;
          Hd.xelem(20) = -r*(1-s2)*(t-1)*t;
          Hd.xelem(21) = -r*(s-1)*s*(1-t2);
          Hd.xelem(22) = ((r+1)*(1-s2)*(1-t2))/2.0E+0+(r*(1-s2)*(1-t2))/2.0E+0;
          Hd.xelem(23) = -r*s*(s+1)*(1-t2);
          Hd.xelem(24) = (r*(1-s2)*(1-t2))/2.0E+0+((r-1)*(1-s2)*(1-t2))/2.0E+0;
          Hd.xelem(25) = -r*(1-s2)*t*(t+1);
          Hd.xelem(26) = -2*r*(1-s2)*(1-t2);
          Hd.xelem(27) = ((r-1)*r*s*(t-1)*t)/8.0E+0+((r-1)*r*(s-1)*(t-1)*t)/8.0E+0;
          Hd.xelem(28) = (r*(r+1)*s*(t-1)*t)/8.0E+0+(r*(r+1)*(s-1)*(t-1)*t)/8.0E+0;
          Hd.xelem(29) = (r*(r+1)*(s+1)*(t-1)*t)/8.0E+0+(r*(r+1)*s*(t-1)*t)/8.0E+0;
          Hd.xelem(30) = ((r-1)*r*(s+1)*(t-1)*t)/8.0E+0+((r-1)*r*s*(t-1)*t)/8.0E+0;
          Hd.xelem(31) = ((r-1)*r*s*t*(t+1))/8.0E+0+((r-1)*r*(s-1)*t*(t+1))/8.0E+0;
          Hd.xelem(32) = (r*(r+1)*s*t*(t+1))/8.0E+0+(r*(r+1)*(s-1)*t*(t+1))/8.0E+0;
          Hd.xelem(33) = (r*(r+1)*(s+1)*t*(t+1))/8.0E+0+(r*(r+1)*s*t*(t+1))/8.0E+0;
          Hd.xelem(34) = ((r-1)*r*(s+1)*t*(t+1))/8.0E+0+((r-1)*r*s*t*(t+1))/8.0E+0;
          Hd.xelem(35) = ((1-r2)*s*(t-1)*t)/4.0E+0+((1-r2)*(s-1)*(t-1)*t)/4.0E+0;
          Hd.xelem(36) = -(r*(r+1)*s*(t-1)*t)/2.0E+0;
          Hd.xelem(37) = ((1-r2)*(s+1)*(t-1)*t)/4.0E+0+((1-r2)*s*(t-1)*t)/4.0E+0;
          Hd.xelem(38) = -((r-1)*r*s*(t-1)*t)/2.0E+0;
          Hd.xelem(39) = ((r-1)*r*s*(1-t2))/4.0E+0+((r-1)*r*(s-1)*(1-t2))/4.0E+0;
          Hd.xelem(40) = (r*(r+1)*s*(1-t2))/4.0E+0+(r*(r+1)*(s-1)*(1-t2))/4.0E+0;
          Hd.xelem(41) = (r*(r+1)*(s+1)*(1-t2))/4.0E+0+(r*(r+1)*s*(1-t2))/4.0E+0;
          Hd.xelem(42) = ((r-1)*r*(s+1)*(1-t2))/4.0E+0+((r-1)*r*s*(1-t2))/4.0E+0;
          Hd.xelem(43) = ((1-r2)*s*t*(t+1))/4.0E+0+((1-r2)*(s-1)*t*(t+1))/4.0E+0;
          Hd.xelem(44) = -(r*(r+1)*s*t*(t+1))/2.0E+0;
          Hd.xelem(45) = ((1-r2)*(s+1)*t*(t+1))/4.0E+0+((1-r2)*s*t*(t+1))/4.0E+0;
          Hd.xelem(46) = -((r-1)*r*s*t*(t+1))/2.0E+0;
          Hd.xelem(47) = -(1-r2)*s*(t-1)*t;
          Hd.xelem(48) = ((1-r2)*s*(1-t2))/2.0E+0+((1-r2)*(s-1)*(1-t2))/2.0E+0;
          Hd.xelem(49) = -r*(r+1)*s*(1-t2);
          Hd.xelem(50) = ((1-r2)*(s+1)*(1-t2))/2.0E+0+((1-r2)*s*(1-t2))/2.0E+0;
          Hd.xelem(51) = -(r-1)*r*s*(1-t2);
          Hd.xelem(52) = -(1-r2)*s*t*(t+1);
          Hd.xelem(53) = -2*(1-r2)*s*(1-t2);
          Hd.xelem(54) = ((r-1)*r*(s-1)*s*t)/8.0E+0+((r-1)*r*(s-1)*s*(t-1))/8.0E+0;
          Hd.xelem(55) = (r*(r+1)*(s-1)*s*t)/8.0E+0+(r*(r+1)*(s-1)*s*(t-1))/8.0E+0;
          Hd.xelem(56) = (r*(r+1)*s*(s+1)*t)/8.0E+0+(r*(r+1)*s*(s+1)*(t-1))/8.0E+0;
          Hd.xelem(57) = ((r-1)*r*s*(s+1)*t)/8.0E+0+((r-1)*r*s*(s+1)*(t-1))/8.0E+0;
          Hd.xelem(58) = ((r-1)*r*(s-1)*s*(t+1))/8.0E+0+((r-1)*r*(s-1)*s*t)/8.0E+0;
          Hd.xelem(59) = (r*(r+1)*(s-1)*s*(t+1))/8.0E+0+(r*(r+1)*(s-1)*s*t)/8.0E+0;
          Hd.xelem(60) = (r*(r+1)*s*(s+1)*(t+1))/8.0E+0+(r*(r+1)*s*(s+1)*t)/8.0E+0;
          Hd.xelem(61) = ((r-1)*r*s*(s+1)*(t+1))/8.0E+0+((r-1)*r*s*(s+1)*t)/8.0E+0;
          Hd.xelem(62) = ((1-r2)*(s-1)*s*t)/4.0E+0+((1-r2)*(s-1)*s*(t-1))/4.0E+0;
          Hd.xelem(63) = (r*(r+1)*(1-s2)*t)/4.0E+0+(r*(r+1)*(1-s2)*(t-1))/4.0E+0;
          Hd.xelem(64) = ((1-r2)*s*(s+1)*t)/4.0E+0+((1-r2)*s*(s+1)*(t-1))/4.0E+0;
          Hd.xelem(65) = ((r-1)*r*(1-s2)*t)/4.0E+0+((r-1)*r*(1-s2)*(t-1))/4.0E+0;
          Hd.xelem(66) = -((r-1)*r*(s-1)*s*t)/2.0E+0;
          Hd.xelem(67) = -(r*(r+1)*(s-1)*s*t)/2.0E+0;
          Hd.xelem(68) = -(r*(r+1)*s*(s+1)*t)/2.0E+0;
          Hd.xelem(69) = -((r-1)*r*s*(s+1)*t)/2.0E+0;
          Hd.xelem(70) = ((1-r2)*(s-1)*s*(t+1))/4.0E+0+((1-r2)*(s-1)*s*t)/4.0E+0;
          Hd.xelem(71) = (r*(r+1)*(1-s2)*(t+1))/4.0E+0+(r*(r+1)*(1-s2)*t)/4.0E+0;
          Hd.xelem(72) = ((1-r2)*s*(s+1)*(t+1))/4.0E+0+((1-r2)*s*(s+1)*t)/4.0E+0;
          Hd.xelem(73) = ((r-1)*r*(1-s2)*(t+1))/4.0E+0+((r-1)*r*(1-s2)*t)/4.0E+0;
          Hd.xelem(74) = ((1-r2)*(1-s2)*t)/2.0E+0+((1-r2)*(1-s2)*(t-1))/2.0E+0;
          Hd.xelem(75) = -(1-r2)*(s-1)*s*t;
          Hd.xelem(76) = -r*(r+1)*(1-s2)*t;
          Hd.xelem(77) = -(1-r2)*s*(s+1)*t;
          Hd.xelem(78) = -(r-1)*r*(1-s2)*t;
          Hd.xelem(79) = ((1-r2)*(1-s2)*(t+1))/2.0E+0+((1-r2)*(1-s2)*t)/2.0E+0;
          Hd.xelem(80) = -2*(1-r2)*(1-s2)*t;
     }

     virtual Matrix InterpGaussToNodal(FemMatrixType eMatType, const Matrix& taug) const override final {
          return InterpGaussToNodalTpl<double>(eMatType, taug);
     }

     virtual ComplexMatrix InterpGaussToNodal(FemMatrixType eMatType, const ComplexMatrix& taug) const override final {
          return InterpGaussToNodalTpl<std::complex<double> >(eMatType, taug);
     }

private:
     template <typename T>
     typename PostProcTypeTraits<T>::MatrixType
     InterpGaussToNodalTpl(FemMatrixType eMatType,
                           const typename PostProcTypeTraits<T>::MatrixType& taug) const {
          const IntegrationRule& oIntegRule = GetIntegrationRule(eMatType);
          const octave_idx_type iNumGauss = oIntegRule.iGetNumEvalPoints();
          const octave_idx_type iNumDir = oIntegRule.iGetNumDirections();
          const octave_idx_type iNumNodes = nodes.numel();
          ColumnVector rv(iNumDir);

          Matrix H(iNumGauss, iNumNodes);

          for (octave_idx_type i = 0; i < iNumGauss; ++i) {
               for (octave_idx_type j = 0; j < iNumDir; ++j) {
                    rv.xelem(j) = oIntegRule.dGetPosition(i, j);
               }

               ScalarInterpMatrix(rv, H, i);
          }

          return H.solve(taug);
     }

     void ScalarInterpMatrix(const ColumnVector& rv, Matrix& Hs, octave_idx_type irow) const override final {
          FEM_ASSERT(rv.numel() == 3);
          FEM_ASSERT(Hs.columns() == 27);
          FEM_ASSERT(irow >= 0);
          FEM_ASSERT(irow < Hs.rows());

          const double r = rv.xelem(0);
          const double s = rv.xelem(1);
          const double t = rv.xelem(2);
          const double r2 = r * r;
          const double s2 = s * s;
          const double t2 = t * t;
          const octave_idx_type nrows = Hs.rows();

          Hs.xelem(irow) = ((r-1)*r*(s-1)*s*(t-1)*t)/8.0E+0;
          Hs.xelem(nrows+irow) = (r*(r+1)*(s-1)*s*(t-1)*t)/8.0E+0;
          Hs.xelem(2*nrows+irow) = (r*(r+1)*s*(s+1)*(t-1)*t)/8.0E+0;
          Hs.xelem(3*nrows+irow) = ((r-1)*r*s*(s+1)*(t-1)*t)/8.0E+0;
          Hs.xelem(4*nrows+irow) = ((r-1)*r*(s-1)*s*t*(t+1))/8.0E+0;
          Hs.xelem(5*nrows+irow) = (r*(r+1)*(s-1)*s*t*(t+1))/8.0E+0;
          Hs.xelem(6*nrows+irow) = (r*(r+1)*s*(s+1)*t*(t+1))/8.0E+0;
          Hs.xelem(7*nrows+irow) = ((r-1)*r*s*(s+1)*t*(t+1))/8.0E+0;
          Hs.xelem(8*nrows+irow) = ((1-r2)*(s-1)*s*(t-1)*t)/4.0E+0;
          Hs.xelem(9*nrows+irow) = (r*(r+1)*(1-s2)*(t-1)*t)/4.0E+0;
          Hs.xelem(10*nrows+irow) = ((1-r2)*s*(s+1)*(t-1)*t)/4.0E+0;
          Hs.xelem(11*nrows+irow) = ((r-1)*r*(1-s2)*(t-1)*t)/4.0E+0;
          Hs.xelem(12*nrows+irow) = ((r-1)*r*(s-1)*s*(1-t2))/4.0E+0;
          Hs.xelem(13*nrows+irow) = (r*(r+1)*(s-1)*s*(1-t2))/4.0E+0;
          Hs.xelem(14*nrows+irow) = (r*(r+1)*s*(s+1)*(1-t2))/4.0E+0;
          Hs.xelem(15*nrows+irow) = ((r-1)*r*s*(s+1)*(1-t2))/4.0E+0;
          Hs.xelem(16*nrows+irow) = ((1-r2)*(s-1)*s*t*(t+1))/4.0E+0;
          Hs.xelem(17*nrows+irow) = (r*(r+1)*(1-s2)*t*(t+1))/4.0E+0;
          Hs.xelem(18*nrows+irow) = ((1-r2)*s*(s+1)*t*(t+1))/4.0E+0;
          Hs.xelem(19*nrows+irow) = ((r-1)*r*(1-s2)*t*(t+1))/4.0E+0;
          Hs.xelem(20*nrows+irow) = ((1-r2)*(1-s2)*(t-1)*t)/2.0E+0;
          Hs.xelem(21*nrows+irow) = ((1-r2)*(s-1)*s*(1-t2))/2.0E+0;
          Hs.xelem(22*nrows+irow) = (r*(r+1)*(1-s2)*(1-t2))/2.0E+0;
          Hs.xelem(23*nrows+irow) = ((1-r2)*s*(s+1)*(1-t2))/2.0E+0;
          Hs.xelem(24*nrows+irow) = ((r-1)*r*(1-s2)*(1-t2))/2.0E+0;
          Hs.xelem(25*nrows+irow) = ((1-r2)*(1-s2)*t*(t+1))/2.0E+0;
          Hs.xelem(26*nrows+irow) = (1-r2)*(1-s2)*(1-t2);
     }

     static octave_idx_type SelectIntegrationRule(FemMatrixType eMatType) {
          switch (eMatType) {
          case MAT_MASS_LUMPED:
               return 1;
          default:
               return 0;
          }
     }

     static array<IntegrationRule, 2> rgIntegRule;
};

array<IntegrationRule, 2> Iso27::rgIntegRule;

class Penta15: public Element3D
{
public:
     Penta15(ElementTypes::TypeId eltype, octave_idx_type id, const Matrix& X, const Material* material, const int32NDArray& nodes, const ElementData& data)
          :Element3D(eltype, id, X, material, nodes, data) {
          FEM_ASSERT(nodes.numel() == 15);
     }

     static void AllocIntegrationRule(FemMatrixType eMatType) {
          // r = 0 ... 1
          // s = 0 ... (1 - r)
          // t = -1 ... 1

          constexpr octave_idx_type N = 3;
          constexpr octave_idx_type M = 7;
          constexpr double r1 = 0.1012865073235;
          constexpr double r2 = 0.7974269853531;
          constexpr double r3 = r1;
          constexpr double r4 = 0.4701420641051;
          constexpr double r5 = r4;
          constexpr double r6 = 0.0597158717898;
          constexpr double r7 = 0.3333333333333;
          constexpr double s1 = r1;
          constexpr double s2 = r1;
          constexpr double s3 = r2;
          constexpr double s4 = r6;
          constexpr double s5 = r4;
          constexpr double s6 = r4;
          constexpr double s7 = r7;
          constexpr double w1 = 0.1259391805448;
          constexpr double w2 = w1;
          constexpr double w3 = w1;
          constexpr double w4 = 0.1323941527885;
          constexpr double w5 = w4;
          constexpr double w6 = w4;
          constexpr double w7 = 0.2250000000001;

          static constexpr double r[2][M] = {{r1, r2, r3, r4, r5, r6, r7}, {0, 1, 0, 0.5,  0.5,    0}};
          static constexpr double s[2][M] = {{s1, s2, s3, s4, s5, s6, s7}, {0, 0, 1,   0,  0.5,  0.5}};
          static constexpr double w[2][M] = {{w1, w2, w3, w4, w5, w6, w7}, {1./6., 1./6., 1./6., 1./6., 1./6., 1./6.}};
          static constexpr double t[2][N] = {{0.774596669241483, 0., -0.774596669241483}, {1., 0., -1.}};
          static constexpr double alpha[2][N] = {{0.555555555555556, 0.888888888888889, 0.555555555555556}, {2./3., 2./3., 2./3.}};

          const IntegRuleType oIRT = SelectIntegrationRule(eMatType);

          if (!rgIntegRule[oIRT.iIntegRule].iGetNumEvalPoints()) {
               rgIntegRule[oIRT.iIntegRule].SetNumEvalPoints(oIRT.iNumPoints * N, 3);

               for (octave_idx_type i = 0; i < oIRT.iNumPoints; ++i) {
                    for (octave_idx_type j = 0; j < N; ++j) {
                         rgIntegRule[oIRT.iIntegRule].SetWeight(i * N + j, 0.5 * w[oIRT.iIntegRule][i] * alpha[oIRT.iIntegRule][j]);
                         rgIntegRule[oIRT.iIntegRule].SetPosition(i * N + j, 0, r[oIRT.iIntegRule][i]);
                         rgIntegRule[oIRT.iIntegRule].SetPosition(i * N + j, 1, s[oIRT.iIntegRule][i]);
                         rgIntegRule[oIRT.iIntegRule].SetPosition(i * N + j, 2, t[oIRT.iIntegRule][j]);
                    }
               }
          }
     }

     virtual const IntegrationRule& GetIntegrationRule(FemMatrixType eMatType) const override final {
          const IntegRuleType oIRT = SelectIntegrationRule(eMatType);

          FEM_ASSERT(oIRT.iIntegRule >= 0);
          FEM_ASSERT(static_cast<size_t>(oIRT.iIntegRule) < rgIntegRule.size());
          FEM_ASSERT(rgIntegRule[oIRT.iIntegRule].iGetNumEvalPoints() > 0);

          return rgIntegRule[oIRT.iIntegRule];
     }

protected:
     virtual void ScalarInterpMatrixDer(const ColumnVector& rv, Matrix& Hd) const override final {
          const double r = rv.xelem(0);
          const double s = rv.xelem(1);
          const double t = rv.xelem(2);
          const double t2 = t * t;

          FEM_ASSERT(rv.numel() == 3);
          FEM_ASSERT(Hd.rows() == 15);
          FEM_ASSERT(Hd.columns() == 3);

          Hd.xelem(0) = (-((1-t)*((-t)-2*s-2*r))/2.0E+0)-((-s)-r+1)*(1-t);
          Hd.xelem(1) = ((1-t)*((-t)+2*r-2))/2.0E+0+r*(1-t);
          Hd.xelem(2) = 0;
          Hd.xelem(3) = (-((t+1)*(t-2*s-2*r))/2.0E+0)-((-s)-r+1)*(t+1);
          Hd.xelem(4) = ((t+1)*(t+2*r-2))/2.0E+0+r*(t+1);
          Hd.xelem(5) = 0;
          Hd.xelem(6) = 2*((-s)-r+1)*(1-t)-2*r*(1-t);
          Hd.xelem(7) = 2*s*(1-t);
          Hd.xelem(8) = -2*s*(1-t);
          Hd.xelem(9) = 2*((-s)-r+1)*(t+1)-2*r*(t+1);
          Hd.xelem(10) = 2*s*(t+1);
          Hd.xelem(11) = -2*s*(t+1);
          Hd.xelem(12) = t2-1;
          Hd.xelem(13) = 1-t2;
          Hd.xelem(14) = 0;
          Hd.xelem(15) = (-((1-t)*((-t)-2*s-2*r))/2.0E+0)-((-s)-r+1)*(1-t);
          Hd.xelem(16) = 0;
          Hd.xelem(17) = ((1-t)*((-t)+2*s-2))/2.0E+0+s*(1-t);
          Hd.xelem(18) = (-((t+1)*(t-2*s-2*r))/2.0E+0)-((-s)-r+1)*(t+1);
          Hd.xelem(19) = 0;
          Hd.xelem(20) = ((t+1)*(t+2*s-2))/2.0E+0+s*(t+1);
          Hd.xelem(21) = -2*r*(1-t);
          Hd.xelem(22) = 2*r*(1-t);
          Hd.xelem(23) = 2*((-s)-r+1)*(1-t)-2*s*(1-t);
          Hd.xelem(24) = -2*r*(t+1);
          Hd.xelem(25) = 2*r*(t+1);
          Hd.xelem(26) = 2*((-s)-r+1)*(t+1)-2*s*(t+1);
          Hd.xelem(27) = t2-1;
          Hd.xelem(28) = 0;
          Hd.xelem(29) = 1-t2;
          Hd.xelem(30) = (-(((-s)-r+1)*((-t)-2*s-2*r))/2.0E+0)-(((-s)-r+1)*(1-t))/2.0E+0;
          Hd.xelem(31) = (-(r*((-t)+2*r-2))/2.0E+0)-(r*(1-t))/2.0E+0;
          Hd.xelem(32) = (-(s*((-t)+2*s-2))/2.0E+0)-(s*(1-t))/2.0E+0;
          Hd.xelem(33) = (((-s)-r+1)*(t-2*s-2*r))/2.0E+0+(((-s)-r+1)*(t+1))/2.0E+0;
          Hd.xelem(34) = (r*(t+2*r-2))/2.0E+0+(r*(t+1))/2.0E+0;
          Hd.xelem(35) = (s*(t+2*s-2))/2.0E+0+(s*(t+1))/2.0E+0;
          Hd.xelem(36) = -2*r*((-s)-r+1);
          Hd.xelem(37) = -2*r*s;
          Hd.xelem(38) = -2*((-s)-r+1)*s;
          Hd.xelem(39) = 2*r*((-s)-r+1);
          Hd.xelem(40) = 2*r*s;
          Hd.xelem(41) = 2*((-s)-r+1)*s;
          Hd.xelem(42) = -2*((-s)-r+1)*t;
          Hd.xelem(43) = -2*r*t;
          Hd.xelem(44) = -2*s*t;
     }

     virtual Matrix InterpGaussToNodal(FemMatrixType eMatType, const Matrix& taug) const override final {
          return InterpGaussToNodalTpl<double>(eMatType, taug);
     }

     virtual ComplexMatrix InterpGaussToNodal(FemMatrixType eMatType, const ComplexMatrix& taug) const override final {
          return InterpGaussToNodalTpl<std::complex<double> >(eMatType, taug);
     }

private:
     template <typename T>
     typename PostProcTypeTraits<T>::MatrixType
     InterpGaussToNodalTpl(FemMatrixType eMatType, const typename PostProcTypeTraits<T>::MatrixType& taug) const {
          const IntegrationRule& oIntegRule = GetIntegrationRule(eMatType);
          const octave_idx_type iNumGauss = oIntegRule.iGetNumEvalPoints();
          const octave_idx_type iNumDir = oIntegRule.iGetNumDirections();
          const octave_idx_type iNumNodes = nodes.numel();
          ColumnVector rv(iNumDir);

          Matrix H(iNumGauss, iNumNodes);

          for (octave_idx_type i = 0; i < iNumGauss; ++i) {
               for (octave_idx_type j = 0; j < iNumDir; ++j) {
                    rv.xelem(j) = oIntegRule.dGetPosition(i, j);
               }

               ScalarInterpMatrix(rv, H, i);
          }

          return H.solve(taug);
     }

     void ScalarInterpMatrix(const ColumnVector& rv, Matrix& Hs, octave_idx_type irow) const override final {
          FEM_ASSERT(rv.numel() == 3);
          FEM_ASSERT(Hs.columns() == 15);
          FEM_ASSERT(irow >= 0);
          FEM_ASSERT(irow < Hs.rows());

          const double r = rv.xelem(0);
          const double s = rv.xelem(1);
          const double t = rv.xelem(2);
          const double t2 = t * t;

          const octave_idx_type nrows = Hs.rows();

          Hs.xelem(irow + nrows * 0) = (((-s)-r+1)*(1-t)*((-t)-2*s-2*r))/2.0E+0;
          Hs.xelem(irow + nrows * 1) = (r*(1-t)*((-t)+2*r-2))/2.0E+0;
          Hs.xelem(irow + nrows * 2) = (s*(1-t)*((-t)+2*s-2))/2.0E+0;
          Hs.xelem(irow + nrows * 3) = (((-s)-r+1)*(t+1)*(t-2*s-2*r))/2.0E+0;
          Hs.xelem(irow + nrows * 4) = (r*(t+1)*(t+2*r-2))/2.0E+0;
          Hs.xelem(irow + nrows * 5) = (s*(t+1)*(t+2*s-2))/2.0E+0;
          Hs.xelem(irow + nrows * 6) = 2*r*((-s)-r+1)*(1-t);
          Hs.xelem(irow + nrows * 7) = 2*r*s*(1-t);
          Hs.xelem(irow + nrows * 8) = 2*((-s)-r+1)*s*(1-t);
          Hs.xelem(irow + nrows * 9) = 2*r*((-s)-r+1)*(t+1);
          Hs.xelem(irow + nrows * 10) = 2*r*s*(t+1);
          Hs.xelem(irow + nrows * 11) = 2*((-s)-r+1)*s*(t+1);
          Hs.xelem(irow + nrows * 12) = ((-s)-r+1)*(1-t2);
          Hs.xelem(irow + nrows * 13) = r*(1-t2);
          Hs.xelem(irow + nrows * 14) = s*(1-t2);
     }

     struct IntegRuleType {
          IntegRuleType(octave_idx_type iIntegRule, octave_idx_type iNumPoints)
               :iIntegRule(iIntegRule), iNumPoints(iNumPoints) {
          }

          octave_idx_type iIntegRule;
          octave_idx_type iNumPoints;
     };

     static IntegRuleType SelectIntegrationRule(FemMatrixType eMatType) {
          switch (eMatType) {
          case MAT_MASS_LUMPED:
               return IntegRuleType{1, 6};
          default:
               return IntegRuleType{0, 7};
          }
     }

     static array<IntegrationRule, 2> rgIntegRule;
};

array<IntegrationRule, 2> Penta15::rgIntegRule;

class Tet10h: public Element3D
{
public:
     Tet10h(ElementTypes::TypeId eltype, octave_idx_type id, const Matrix& X, const Material* material, const int32NDArray& nodes, const ElementData& data)
          :Element3D(eltype, id, X, material, nodes, data) {
          FEM_ASSERT(nodes.numel() == 10);
     }

     static void AllocIntegrationRule(FemMatrixType eMatType) {
          switch (eMatType) {
          case MAT_DAMPING:
          case MAT_DAMPING_SYM:
          case MAT_DAMPING_SYM_L:
          case MAT_DAMPING_FLUID_STRUCT_RE:
          case MAT_DAMPING_FLUID_STRUCT_IM:
               AllocIntegrationRule(MAT_STIFFNESS);
               AllocIntegrationRule(MAT_MASS);
               return;

          default:
               break;
          }

          constexpr double a1 = 0.25;
          constexpr double b1 = 1. / 6.;
          constexpr double c1 = 0.5;
          constexpr double d1 = -2. / 15.;
          constexpr double e1 = 3. / 40.;
          constexpr octave_idx_type N1 = 5;
          static constexpr double r1[N1] = {a1, b1, b1, b1, c1};
          static constexpr double s1[N1] = {a1, b1, b1, c1, b1};
          static constexpr double t1[N1] = {a1, b1, c1, b1, b1};
          static constexpr double w1[N1] = {d1, e1, e1, e1, e1};

          constexpr double a2 = 0.25;
          constexpr double b2_1 = (7. + sqrt(15.)) / 34.;
          constexpr double b2_2 = (7. - sqrt(15.)) / 34.;
          constexpr double c2_1 = (13. - 3. * sqrt(15.)) / 34.;
          constexpr double c2_2 = (13. + 3. * sqrt(15.)) / 34.;
          constexpr double d2 = (5. - sqrt(15.)) / 20.;
          constexpr double e2 = (5. + sqrt(15.)) / 20.;
          constexpr double f2 = 8. / 405.;
          constexpr double g2 = (2665. - 14. * sqrt(15.)) / 226800.;
          constexpr double h2 = (2665. + 14. * sqrt(15.)) / 226800.;
          constexpr double i2 = 5. / 567.;
          constexpr octave_idx_type N2 = 15;

          static constexpr double t2[N2] = {a2, b2_1, b2_1, b2_1, c2_1, b2_2, b2_2, b2_2, c2_2, d2, d2, e2, d2, e2, e2};

          static constexpr double r2[N2] = {a2, b2_1, b2_1, c2_1, b2_1, b2_2, b2_2, c2_2, b2_2, d2, e2, d2, e2, d2, e2};

          static constexpr double s2[N2] = {a2, b2_1, c2_1, b2_1, b2_1, b2_2, c2_2, b2_2, b2_2, e2, d2, d2, e2, e2, d2};

          static constexpr double w2[N2] = {f2, g2, g2, g2, g2, h2, h2, h2, h2, i2, i2, i2, i2, i2, i2};

          constexpr double a3 = (5. - sqrt(5.)) / 20.;
          constexpr double b3 = (5. + 3. * sqrt(5)) / 20.;
          constexpr double c3 = 1. / 24.;
          constexpr octave_idx_type N3 = 4;

          static constexpr double r3[N3] = {a3, a3, a3, b3};
          static constexpr double s3[N3] = {a3, a3, b3, a3};
          static constexpr double t3[N3] = {a3, b3, a3, a3};
          static constexpr double w3[N3] = {c3, c3, c3, c3};

          static constexpr octave_idx_type Ni[RNUM] = {N1, N2, N3};
          static constexpr const double* ri[RNUM] = {r1, r2, r3};
          static constexpr const double* si[RNUM] = {s1, s2, s3};
          static constexpr const double* ti[RNUM] = {t1, t2, t3};
          static constexpr const double* wi[RNUM] = {w1, w2, w3};

          const octave_idx_type iIntegRule = SelectIntegrationRule(eMatType);

          if (!rgIntegRule[iIntegRule].iGetNumEvalPoints()) {
               const octave_idx_type N = Ni[iIntegRule];
               const double* const r = ri[iIntegRule];
               const double* const s = si[iIntegRule];
               const double* const t = ti[iIntegRule];
               const double* const w = wi[iIntegRule];

               rgIntegRule[iIntegRule].SetNumEvalPoints(N, 3);

               for (octave_idx_type i = 0; i < N; ++i) {
                    rgIntegRule[iIntegRule].SetWeight(i, w[i]);
                    rgIntegRule[iIntegRule].SetPosition(i, 0, r[i]);
                    rgIntegRule[iIntegRule].SetPosition(i, 1, s[i]);
                    rgIntegRule[iIntegRule].SetPosition(i, 2, t[i]);
               }
          }
     }

     virtual const IntegrationRule& GetIntegrationRule(FemMatrixType eMatType) const override final {
          const octave_idx_type iIntegRule = SelectIntegrationRule(eMatType);

          FEM_ASSERT(iIntegRule >= 0);
          FEM_ASSERT(static_cast<size_t>(iIntegRule) < rgIntegRule.size());
          FEM_ASSERT(rgIntegRule[iIntegRule].iGetNumEvalPoints() > 0);

          return rgIntegRule[iIntegRule];
     }

protected:
     virtual Matrix InterpGaussToNodal(FemMatrixType eMatType, const Matrix& taug) const override final {
          return InterpGaussToNodalTpl<double>(eMatType, taug);
     }

     virtual ComplexMatrix InterpGaussToNodal(FemMatrixType eMatType, const ComplexMatrix& taug) const override final {
          return InterpGaussToNodalTpl<std::complex<double> >(eMatType, taug);
     }

     void ScalarInterpMatrix(const ColumnVector& rv, Matrix& Hs, octave_idx_type irow) const override final {
          FEM_ASSERT(rv.numel() == 3);
          FEM_ASSERT(Hs.columns() == 10);
          FEM_ASSERT(irow >= 0);
          FEM_ASSERT(irow < Hs.rows());

          const double r = rv.xelem(0);
          const double s = rv.xelem(1);
          const double t = rv.xelem(2);
          const octave_idx_type nrows = Hs.rows();

          Hs.xelem(irow + nrows * 0) = s*(2*s-1);
          Hs.xelem(irow + nrows * 1) = t*(2*t-1);
          Hs.xelem(irow + nrows * 2) = ((-2*t)-2*s-2*r+1)*((-t)-s-r+1);
          Hs.xelem(irow + nrows * 3) = r*(2*r-1);
          Hs.xelem(irow + nrows * 4) = 4*s*t;
          Hs.xelem(irow + nrows * 5) = 4*((-t)-s-r+1)*t;
          Hs.xelem(irow + nrows * 6) = 4*s*((-t)-s-r+1);
          Hs.xelem(irow + nrows * 7) = 4*r*s;
          Hs.xelem(irow + nrows * 8) = 4*r*t;
          Hs.xelem(irow + nrows * 9) = 4*r*((-t)-s-r+1);
     }

     virtual void ScalarInterpMatrixDer(const ColumnVector& rv, Matrix& Hd) const override final {
          const double r = rv.xelem(0);
          const double s = rv.xelem(1);
          const double t = rv.xelem(2);

          FEM_ASSERT(rv.numel() == 3);
          FEM_ASSERT(Hd.rows() == 10);
          FEM_ASSERT(Hd.columns() == 3);

          Hd.xelem(0) = 0;
          Hd.xelem(1) = 0;
          Hd.xelem(2) = 2*t-2*((-t)-s-r+1)+2*s+2*r-1;
          Hd.xelem(3) = 4*r-1;
          Hd.xelem(4) = 0;
          Hd.xelem(5) = -4*t;
          Hd.xelem(6) = -4*s;
          Hd.xelem(7) = 4*s;
          Hd.xelem(8) = 4*t;
          Hd.xelem(9) = 4*((-t)-s-r+1)-4*r;
          Hd.xelem(10) = 4*s-1;
          Hd.xelem(11) = 0;
          Hd.xelem(12) = 2*t-2*((-t)-s-r+1)+2*s+2*r-1;
          Hd.xelem(13) = 0;
          Hd.xelem(14) = 4*t;
          Hd.xelem(15) = -4*t;
          Hd.xelem(16) = 4*((-t)-s-r+1)-4*s;
          Hd.xelem(17) = 4*r;
          Hd.xelem(18) = 0;
          Hd.xelem(19) = -4*r;
          Hd.xelem(20) = 0;
          Hd.xelem(21) = 4*t-1;
          Hd.xelem(22) = 2*t-2*((-t)-s-r+1)+2*s+2*r-1;
          Hd.xelem(23) = 0;
          Hd.xelem(24) = 4*s;
          Hd.xelem(25) = 4*((-t)-s-r+1)-4*t;
          Hd.xelem(26) = -4*s;
          Hd.xelem(27) = 0;
          Hd.xelem(28) = 4*r;
          Hd.xelem(29) = -4*r;
     }
private:
     enum IntegRuleType {
          R1 = 0,
          R2,
          R3,
          RNUM
     };

     template <typename T>
     typename PostProcTypeTraits<T>::MatrixType
     InterpGaussToNodalTpl(FemMatrixType eMatType,
                           const typename PostProcTypeTraits<T>::MatrixType& taug) const {
          const IntegrationRule& oIntegRule = GetIntegrationRule(eMatType);
          const octave_idx_type iNumGauss = oIntegRule.iGetNumEvalPoints();
          const octave_idx_type iNumDir = oIntegRule.iGetNumDirections();
          const octave_idx_type iNumNodesCorner = 4;
          ColumnVector rv(iNumDir);

          Matrix H(iNumGauss, iNumNodesCorner);

          for (octave_idx_type i = 0; i < iNumGauss; ++i) {
               for (octave_idx_type j = 0; j < iNumDir; ++j) {
                    rv.xelem(j) = oIntegRule.dGetPosition(i, j);
               }

               ScalarInterpMatrixCornerNodes(rv, H, i);
          }

          typename PostProcTypeTraits<T>::MatrixType taun = H.solve(taug);

          taun.resize(nodes.numel(), taun.columns());

          static constexpr struct {
               octave_idx_type ico1, ico2, imid;
          } idxint[] = {
               {0, 1, 4},
               {1, 2, 5},
               {2, 0, 6},
               {0, 3, 7},
               {1, 3, 8},
               {2, 3, 9}
          };

          for (octave_idx_type j = 0; j < taun.columns(); ++j) {
               for (const auto& idx:idxint) {
                    taun.xelem(idx.imid, j) = 0.5 * (taun.xelem(idx.ico1, j) + taun.xelem(idx.ico2, j));
               }
          }

          return taun;
     }

     void ScalarInterpMatrixCornerNodes(const ColumnVector& rv, Matrix& Hs, octave_idx_type irow) const {
          FEM_ASSERT(rv.numel() == 3);
          FEM_ASSERT(Hs.columns() == 4);
          FEM_ASSERT(irow >= 0);
          FEM_ASSERT(irow < Hs.rows());

          const double r = rv.xelem(0);
          const double s = rv.xelem(1);
          const double t = rv.xelem(2);

          const octave_idx_type nrows = Hs.rows();

          Hs.xelem(irow + nrows *  0) = s;
          Hs.xelem(irow + nrows *  1) = t;
          Hs.xelem(irow + nrows *  2) = 1. - r - s - t;
          Hs.xelem(irow + nrows *  3) = r;
     }

     static octave_idx_type SelectIntegrationRule(FemMatrixType eMatType) {
          switch (eMatType & ~MAT_COLL_PNT_OUTPUT) {
          case MAT_STIFFNESS:
          case MAT_STIFFNESS_SYM:
          case MAT_STIFFNESS_SYM_L:
          case MAT_STIFFNESS_TAU0:
          case MAT_STIFFNESS_FLUID_STRUCT_RE:
          case MAT_STIFFNESS_FLUID_STRUCT_IM:
          case VEC_LOAD_CONSISTENT:
          case VEC_LOAD_LUMPED:
          case VEC_LOAD_THERMAL:
          case VEC_LOAD_ACOUSTICS:
          case VEC_LOAD_FLUID_STRUCT:
          case MAT_THERMAL_COND:
          case MAT_STIFFNESS_ACOUSTICS_RE:
          case MAT_STIFFNESS_ACOUSTICS_IM:
          case MAT_DAMPING_ACOUSTICS_RE:
          case MAT_DAMPING_ACOUSTICS_IM:
          case MAT_DAMPING_FLUID_STRUCT_RE:
          case MAT_DAMPING_FLUID_STRUCT_IM:
               return R1;

          case MAT_MASS:
          case MAT_MASS_SYM:
          case MAT_MASS_SYM_L:
          case MAT_STIFFNESS_OMEGA:
          case MAT_STIFFNESS_OMEGA_DOT:
          case MAT_DAMPING_OMEGA:
          case MAT_MASS_FLUID_STRUCT_RE:
          case MAT_MASS_FLUID_STRUCT_IM:
          case VEC_INERTIA_M1:
          case MAT_INERTIA_J:
          case MAT_INERTIA_INV3:
          case MAT_INERTIA_INV4:
          case MAT_INERTIA_INV5:
          case MAT_INERTIA_INV8:
          case MAT_INERTIA_INV9:
          case SCA_TOT_MASS:
          case MAT_HEAT_CAPACITY:
          case MAT_MASS_ACOUSTICS_RE:
          case MAT_MASS_ACOUSTICS_IM:
               return R2;

          case VEC_STRESS_CAUCH:
          case VEC_STRAIN_TOTAL:
          case SCA_STRESS_VMIS:
          case VEC_PARTICLE_VELOCITY:
          case VEC_PARTICLE_VELOCITY_C:
               return R3;

          default:
               throw std::runtime_error("tet10h: matrix type not supported");
          }
     }

     static array<IntegrationRule, RNUM> rgIntegRule;
};

array<IntegrationRule, Tet10h::RNUM> Tet10h::rgIntegRule;

class Tet10: public Element3D
{
     static constexpr double gamma = 1. / 6.;

public:
     Tet10(ElementTypes::TypeId eltype, octave_idx_type id, const Matrix& X, const Material* material, const int32NDArray& nodes, const ElementData& data)
          :Element3D(eltype, id, X, material, nodes, data) {
          FEM_ASSERT(nodes.numel() == 10);
     }

     static void AllocIntegrationRule(FemMatrixType eMatType) {
          switch (eMatType & ~MAT_COLL_PNT_OUTPUT) {
          case MAT_STIFFNESS:
          case MAT_STIFFNESS_SYM:
          case MAT_STIFFNESS_SYM_L:
          case MAT_STIFFNESS_TAU0:
          case VEC_STRESS_CAUCH:
          case VEC_STRAIN_TOTAL:
          case VEC_LOAD_CONSISTENT:
          case VEC_LOAD_LUMPED:
          case VEC_LOAD_THERMAL:
          case VEC_LOAD_FLUID_STRUCT:
          case VEC_LOAD_ACOUSTICS:
          case MAT_THERMAL_COND:
          case MAT_STIFFNESS_ACOUSTICS_RE:
          case MAT_STIFFNESS_ACOUSTICS_IM:
          case MAT_DAMPING_ACOUSTICS_RE:
          case MAT_DAMPING_ACOUSTICS_IM:
          case MAT_STIFFNESS_FLUID_STRUCT_RE:
          case MAT_STIFFNESS_FLUID_STRUCT_IM:
          case VEC_PARTICLE_VELOCITY:
          case VEC_PARTICLE_VELOCITY_C:
               if (!oIntegStiff.iGetNumEvalPoints()) {
                    constexpr double alpha = (5. + 3. * sqrt(5.)) / 20.;
                    constexpr double beta = (5. - sqrt(5.)) / 20.;

                    oIntegStiff.SetNumEvalPoints(4, 4);

                    for (octave_idx_type i = 0; i < 4; ++i) {
                         oIntegStiff.SetWeight(i, 0.25);

                         for (octave_idx_type j = 0; j < 4; ++j) {
                              oIntegStiff.SetPosition(i, j, i == j ? alpha : beta);
                         }
                    }
               }
               break;

          case MAT_MASS:
          case MAT_MASS_SYM:
          case MAT_MASS_SYM_L:
          case MAT_STIFFNESS_OMEGA:
          case MAT_STIFFNESS_OMEGA_DOT:
          case MAT_DAMPING_OMEGA:
          case MAT_MASS_FLUID_STRUCT_RE:
          case MAT_MASS_FLUID_STRUCT_IM:
          case VEC_INERTIA_M1:
          case MAT_INERTIA_J:
          case MAT_INERTIA_INV3:
          case MAT_INERTIA_INV4:
          case MAT_INERTIA_INV5:
          case MAT_INERTIA_INV8:
          case MAT_INERTIA_INV9:
          case MAT_HEAT_CAPACITY:
          case MAT_MASS_ACOUSTICS_RE:
          case MAT_MASS_ACOUSTICS_IM:
               if (!oIntegMass.iGetNumEvalPoints()) {
                    constexpr double g1 = 0.09273525031089122640232391373703060;
                    constexpr double g2 = 0.31088591926330060979734573376345783;
                    constexpr double g3 = 0.45449629587435035050811947372066056;
                    constexpr double w1 = (-1+6*g2*(2+g2*(-7+8*g2))+14*g3-60*g2*(3+4*g2*(-3+4*g2))*g3+4*(-7+30*g2*(3+4*g2*(-3+4*g2)))*g3*g3)/(120*(g1-g2)*(g2*(-3+8*g2)+6*g3+8*g2*(-3+4*g2)*g3-4*(3+4*g2*(-3+4*g2))*g3*g3+8*g1*g1*(1+12*g2*(-1+2*g2)+4*g3-8*g3*g3)+g1*(-3-96*g2*g2+24*g3*(-1+2*g3)+g2*(44+32*(1-2*g3)*g3))));
                    constexpr double w2 = (-1-20*(1+12*g1*(2*g1-1))*w1+20*g3*(2*g3-1)*(4*w1-1))/(20*(1+12*g2*(2*g2-1)+4*g3-8*g3*g3));
                    static constexpr octave_idx_type jk6[][2] = {{1 - 1, 2 - 1}, {1 - 1, 3 - 1}, {1 - 1, 4 - 1}, {2 - 1, 3 - 1}, {2 - 1, 4 - 1}, {3 - 1, 4 - 1}};

                    oIntegMass.SetNumEvalPoints(14, 4);

                    for (octave_idx_type i = 0; i < 5 - 1; ++i) {
                         oIntegMass.SetWeight(i, w1);

                         for (octave_idx_type j = 0; j < 4; ++j) {
                              oIntegMass.SetPosition(i, j, j == i ?  1 - 3 * g1 : g1);
                         }
                    }

                    for (octave_idx_type i = 5 - 1; i < 9 - 1; ++i) {
                         oIntegMass.SetWeight(i, w2);

                         for (octave_idx_type j = 0; j < 4; ++j) {
                              oIntegMass.SetPosition(i, j, j == i - 4 ?  1 - 3 * g2 : g2);
                         }
                    }

                    for (octave_idx_type i = 9 - 1; i < 14; ++i) {
                         oIntegMass.SetWeight(i, 1 / 6. - 2. * (w1 + w2) / 3.);

                         for (octave_idx_type j = 0; j < 4; ++j) {
                              oIntegMass.SetPosition(i, j, (j == jk6[i - 8][0] || j == jk6[i - 8][1]) ? 1. / 2. - g3 : g3);
                         }
                    }
               }
               break;

          case MAT_MASS_LUMPED:
               if (!oIntegMassDiag.iGetNumEvalPoints()) {
                    oIntegMassDiag.SetNumEvalPoints(10, 4);

                    constexpr double w1 = 1. / 10.;
                    constexpr double alpha = 1.;
                    constexpr double beta = 0.;

                    for (octave_idx_type i = 0; i < 4; ++i) {
                         oIntegMassDiag.SetWeight(i, w1);

                         for (octave_idx_type j = 0; j < 4; ++j) {
                              oIntegMassDiag.SetPosition(i, j, i == j ? alpha : beta);
                         }
                    }

                    constexpr double w2 = 1. / 10.;
                    static constexpr double g2[][4] = {{0.5,0.5,0,0},
                                                       {0.5,0,0.5,0},
                                                       {0.5,0,0,0.5},
                                                       {0,0.5,0.5,0},
                                                       {0,0.5,0,0.5},
                                                       {0,0,0.5,0.5}};

                    for (octave_idx_type i = 0; i < 6; ++i) {
                         oIntegMassDiag.SetWeight(i + 4, w2);

                         for (octave_idx_type j = 0; j < 4; ++j) {
                              oIntegMassDiag.SetPosition(i + 4, j, g2[i][j]);
                         }
                    }
               }
               break;

          case MAT_DAMPING:
          case MAT_DAMPING_SYM:
          case MAT_DAMPING_SYM_L:
          case MAT_DAMPING_FLUID_STRUCT_RE:
          case MAT_DAMPING_FLUID_STRUCT_IM:
               AllocIntegrationRule(MAT_MASS);
               AllocIntegrationRule(MAT_STIFFNESS);
               break;

          default:
               throw std::runtime_error("tet10: unknown matrix type");
          }
     }

     virtual const IntegrationRule& GetIntegrationRule(FemMatrixType eMatType) const override final {
          const IntegrationRule& oIntegRule = SelectIntegrationRule(eMatType);

          FEM_ASSERT(oIntegRule.iGetNumEvalPoints() > 0);

          return oIntegRule;
     }

protected:
     virtual double Jacobian(const ColumnVector& rv, Matrix& J) const override final {
          FEM_ASSERT(J.rows() == 4);
          FEM_ASSERT(J.columns() == 4);
          FEM_ASSERT(rv.numel() == 4);

          const double Zeta1 = rv.xelem(0);
          const double Zeta2 = rv.xelem(1);
          const double Zeta3 = rv.xelem(2);
          const double Zeta4 = rv.xelem(3);

          FEM_ASSERT(Zeta1 + Zeta2 + Zeta3 + Zeta4 == 1);

          J.xelem(0) = 1;
          J.xelem(1) = 4*X.xelem(21)*Zeta4+4*X.xelem(18)*Zeta3+4*X.xelem(12)*Zeta2+X.xelem(0)*(4*Zeta1-1);
          J.xelem(2) = 4*X.xelem(22)*Zeta4+4*X.xelem(19)*Zeta3+4*X.xelem(13)*Zeta2+X.xelem(1)*(4*Zeta1-1);
          J.xelem(3) = 4*X.xelem(23)*Zeta4+4*X.xelem(20)*Zeta3+4*X.xelem(14)*Zeta2+X.xelem(2)*(4*Zeta1-1);
          J.xelem(4) = 1;
          J.xelem(5) = 4*X.xelem(24)*Zeta4+4*X.xelem(15)*Zeta3+X.xelem(3)*(4*Zeta2-1)+4*X.xelem(12)*Zeta1;
          J.xelem(6) = 4*X.xelem(25)*Zeta4+4*X.xelem(16)*Zeta3+X.xelem(4)*(4*Zeta2-1)+4*X.xelem(13)*Zeta1;
          J.xelem(7) = 4*X.xelem(26)*Zeta4+4*X.xelem(17)*Zeta3+X.xelem(5)*(4*Zeta2-1)+4*X.xelem(14)*Zeta1;
          J.xelem(8) = 1;
          J.xelem(9) = 4*X.xelem(27)*Zeta4+X.xelem(6)*(4*Zeta3-1)+4*X.xelem(15)*Zeta2+4*X.xelem(18)*Zeta1;
          J.xelem(10) = 4*X.xelem(28)*Zeta4+X.xelem(7)*(4*Zeta3-1)+4*X.xelem(16)*Zeta2+4*X.xelem(19)*Zeta1;
          J.xelem(11) = 4*X.xelem(29)*Zeta4+X.xelem(8)*(4*Zeta3-1)+4*X.xelem(17)*Zeta2+4*X.xelem(20)*Zeta1;
          J.xelem(12) = 1;
          J.xelem(13) = X.xelem(9)*(4*Zeta4-1)+4*X.xelem(27)*Zeta3+4*X.xelem(24)*Zeta2+4*X.xelem(21)*Zeta1;
          J.xelem(14) = X.xelem(10)*(4*Zeta4-1)+4*X.xelem(28)*Zeta3+4*X.xelem(25)*Zeta2+4*X.xelem(22)*Zeta1;
          J.xelem(15) = X.xelem(11)*(4*Zeta4-1)+4*X.xelem(29)*Zeta3+4*X.xelem(26)*Zeta2+4*X.xelem(23)*Zeta1;

          return Determinant4x4(J) * gamma;
     }

     virtual void ScalarGradientMatrix(const ColumnVector& rv, const Matrix& J, const double detJ, Matrix& invJ, Matrix& Bt) const override final {
          FEM_ASSERT(Bt.rows() == 3);
          FEM_ASSERT(Bt.columns() == 10);
          FEM_ASSERT(rv.numel() == 4);
          FEM_ASSERT(J.rows() == 4);
          FEM_ASSERT(J.columns() == 4);
          FEM_ASSERT(invJ.rows() == 4);
          FEM_ASSERT(invJ.columns() == 4);

          const double Zeta1 = rv.xelem(0);
          const double Zeta2 = rv.xelem(1);
          const double Zeta3 = rv.xelem(2);
          const double Zeta4 = rv.xelem(3);

          Inverse4x4(J, detJ / gamma, invJ);

          Bt.xelem(0) = invJ.xelem(4)*(4*Zeta1-1);
          Bt.xelem(1) = invJ.xelem(8)*(4*Zeta1-1);
          Bt.xelem(2) = invJ.xelem(12)*(4*Zeta1-1);
          Bt.xelem(3) = invJ.xelem(5)*(4*Zeta2-1);
          Bt.xelem(4) = invJ.xelem(9)*(4*Zeta2-1);
          Bt.xelem(5) = invJ.xelem(13)*(4*Zeta2-1);
          Bt.xelem(6) = invJ.xelem(6)*(4*Zeta3-1);
          Bt.xelem(7) = invJ.xelem(10)*(4*Zeta3-1);
          Bt.xelem(8) = invJ.xelem(14)*(4*Zeta3-1);
          Bt.xelem(9) = invJ.xelem(7)*(4*Zeta4-1);
          Bt.xelem(10) = invJ.xelem(11)*(4*Zeta4-1);
          Bt.xelem(11) = invJ.xelem(15)*(4*Zeta4-1);
          Bt.xelem(12) = 4*invJ.xelem(4)*Zeta2+4*invJ.xelem(5)*Zeta1;
          Bt.xelem(13) = 4*invJ.xelem(8)*Zeta2+4*invJ.xelem(9)*Zeta1;
          Bt.xelem(14) = 4*invJ.xelem(12)*Zeta2+4*invJ.xelem(13)*Zeta1;
          Bt.xelem(15) = 4*invJ.xelem(5)*Zeta3+4*invJ.xelem(6)*Zeta2;
          Bt.xelem(16) = 4*invJ.xelem(9)*Zeta3+4*invJ.xelem(10)*Zeta2;
          Bt.xelem(17) = 4*invJ.xelem(13)*Zeta3+4*invJ.xelem(14)*Zeta2;
          Bt.xelem(18) = 4*invJ.xelem(4)*Zeta3+4*invJ.xelem(6)*Zeta1;
          Bt.xelem(19) = 4*invJ.xelem(8)*Zeta3+4*invJ.xelem(10)*Zeta1;
          Bt.xelem(20) = 4*invJ.xelem(12)*Zeta3+4*invJ.xelem(14)*Zeta1;
          Bt.xelem(21) = 4*invJ.xelem(4)*Zeta4+4*invJ.xelem(7)*Zeta1;
          Bt.xelem(22) = 4*invJ.xelem(8)*Zeta4+4*invJ.xelem(11)*Zeta1;
          Bt.xelem(23) = 4*invJ.xelem(12)*Zeta4+4*invJ.xelem(15)*Zeta1;
          Bt.xelem(24) = 4*invJ.xelem(5)*Zeta4+4*invJ.xelem(7)*Zeta2;
          Bt.xelem(25) = 4*invJ.xelem(9)*Zeta4+4*invJ.xelem(11)*Zeta2;
          Bt.xelem(26) = 4*invJ.xelem(13)*Zeta4+4*invJ.xelem(15)*Zeta2;
          Bt.xelem(27) = 4*invJ.xelem(6)*Zeta4+4*invJ.xelem(7)*Zeta3;
          Bt.xelem(28) = 4*invJ.xelem(10)*Zeta4+4*invJ.xelem(11)*Zeta3;
          Bt.xelem(29) = 4*invJ.xelem(14)*Zeta4+4*invJ.xelem(15)*Zeta3;
     }

     virtual void StrainMatrix(const ColumnVector& rv, const Matrix& J, const double detJ, Matrix& invJ, Matrix& B) const override final {
          FEM_ASSERT(J.rows() == 4);
          FEM_ASSERT(J.columns() == 4);
          FEM_ASSERT(invJ.rows() == 4);
          FEM_ASSERT(invJ.columns() == 4);
          FEM_ASSERT(B.rows() == 6);
          FEM_ASSERT(B.columns() == 30);
          FEM_ASSERT(rv.numel() == 4);

          const double Zeta1 = rv.xelem(0);
          const double Zeta2 = rv.xelem(1);
          const double Zeta3 = rv.xelem(2);
          const double Zeta4 = rv.xelem(3);

          Inverse4x4(J, detJ / gamma, invJ);

          B.xelem(0) = invJ.xelem(4)*(4*Zeta1-1);
          B.xelem(1) = 0;
          B.xelem(2) = 0;
          B.xelem(3) = invJ.xelem(8)*(4*Zeta1-1);
          B.xelem(4) = 0;
          B.xelem(5) = invJ.xelem(12)*(4*Zeta1-1);
          B.xelem(6) = 0;
          B.xelem(7) = invJ.xelem(8)*(4*Zeta1-1);
          B.xelem(8) = 0;
          B.xelem(9) = invJ.xelem(4)*(4*Zeta1-1);
          B.xelem(10) = invJ.xelem(12)*(4*Zeta1-1);
          B.xelem(11) = 0;
          B.xelem(12) = 0;
          B.xelem(13) = 0;
          B.xelem(14) = invJ.xelem(12)*(4*Zeta1-1);
          B.xelem(15) = 0;
          B.xelem(16) = invJ.xelem(8)*(4*Zeta1-1);
          B.xelem(17) = invJ.xelem(4)*(4*Zeta1-1);
          B.xelem(18) = invJ.xelem(5)*(4*Zeta2-1);
          B.xelem(19) = 0;
          B.xelem(20) = 0;
          B.xelem(21) = invJ.xelem(9)*(4*Zeta2-1);
          B.xelem(22) = 0;
          B.xelem(23) = invJ.xelem(13)*(4*Zeta2-1);
          B.xelem(24) = 0;
          B.xelem(25) = invJ.xelem(9)*(4*Zeta2-1);
          B.xelem(26) = 0;
          B.xelem(27) = invJ.xelem(5)*(4*Zeta2-1);
          B.xelem(28) = invJ.xelem(13)*(4*Zeta2-1);
          B.xelem(29) = 0;
          B.xelem(30) = 0;
          B.xelem(31) = 0;
          B.xelem(32) = invJ.xelem(13)*(4*Zeta2-1);
          B.xelem(33) = 0;
          B.xelem(34) = invJ.xelem(9)*(4*Zeta2-1);
          B.xelem(35) = invJ.xelem(5)*(4*Zeta2-1);
          B.xelem(36) = invJ.xelem(6)*(4*Zeta3-1);
          B.xelem(37) = 0;
          B.xelem(38) = 0;
          B.xelem(39) = invJ.xelem(10)*(4*Zeta3-1);
          B.xelem(40) = 0;
          B.xelem(41) = invJ.xelem(14)*(4*Zeta3-1);
          B.xelem(42) = 0;
          B.xelem(43) = invJ.xelem(10)*(4*Zeta3-1);
          B.xelem(44) = 0;
          B.xelem(45) = invJ.xelem(6)*(4*Zeta3-1);
          B.xelem(46) = invJ.xelem(14)*(4*Zeta3-1);
          B.xelem(47) = 0;
          B.xelem(48) = 0;
          B.xelem(49) = 0;
          B.xelem(50) = invJ.xelem(14)*(4*Zeta3-1);
          B.xelem(51) = 0;
          B.xelem(52) = invJ.xelem(10)*(4*Zeta3-1);
          B.xelem(53) = invJ.xelem(6)*(4*Zeta3-1);
          B.xelem(54) = invJ.xelem(7)*(4*Zeta4-1);
          B.xelem(55) = 0;
          B.xelem(56) = 0;
          B.xelem(57) = invJ.xelem(11)*(4*Zeta4-1);
          B.xelem(58) = 0;
          B.xelem(59) = invJ.xelem(15)*(4*Zeta4-1);
          B.xelem(60) = 0;
          B.xelem(61) = invJ.xelem(11)*(4*Zeta4-1);
          B.xelem(62) = 0;
          B.xelem(63) = invJ.xelem(7)*(4*Zeta4-1);
          B.xelem(64) = invJ.xelem(15)*(4*Zeta4-1);
          B.xelem(65) = 0;
          B.xelem(66) = 0;
          B.xelem(67) = 0;
          B.xelem(68) = invJ.xelem(15)*(4*Zeta4-1);
          B.xelem(69) = 0;
          B.xelem(70) = invJ.xelem(11)*(4*Zeta4-1);
          B.xelem(71) = invJ.xelem(7)*(4*Zeta4-1);
          B.xelem(72) = 4*invJ.xelem(4)*Zeta2+4*invJ.xelem(5)*Zeta1;
          B.xelem(73) = 0;
          B.xelem(74) = 0;
          B.xelem(75) = 4*invJ.xelem(8)*Zeta2+4*invJ.xelem(9)*Zeta1;
          B.xelem(76) = 0;
          B.xelem(77) = 4*invJ.xelem(12)*Zeta2+4*invJ.xelem(13)*Zeta1;
          B.xelem(78) = 0;
          B.xelem(79) = 4*invJ.xelem(8)*Zeta2+4*invJ.xelem(9)*Zeta1;
          B.xelem(80) = 0;
          B.xelem(81) = 4*invJ.xelem(4)*Zeta2+4*invJ.xelem(5)*Zeta1;
          B.xelem(82) = 4*invJ.xelem(12)*Zeta2+4*invJ.xelem(13)*Zeta1;
          B.xelem(83) = 0;
          B.xelem(84) = 0;
          B.xelem(85) = 0;
          B.xelem(86) = 4*invJ.xelem(12)*Zeta2+4*invJ.xelem(13)*Zeta1;
          B.xelem(87) = 0;
          B.xelem(88) = 4*invJ.xelem(8)*Zeta2+4*invJ.xelem(9)*Zeta1;
          B.xelem(89) = 4*invJ.xelem(4)*Zeta2+4*invJ.xelem(5)*Zeta1;
          B.xelem(90) = 4*invJ.xelem(5)*Zeta3+4*invJ.xelem(6)*Zeta2;
          B.xelem(91) = 0;
          B.xelem(92) = 0;
          B.xelem(93) = 4*invJ.xelem(9)*Zeta3+4*invJ.xelem(10)*Zeta2;
          B.xelem(94) = 0;
          B.xelem(95) = 4*invJ.xelem(13)*Zeta3+4*invJ.xelem(14)*Zeta2;
          B.xelem(96) = 0;
          B.xelem(97) = 4*invJ.xelem(9)*Zeta3+4*invJ.xelem(10)*Zeta2;
          B.xelem(98) = 0;
          B.xelem(99) = 4*invJ.xelem(5)*Zeta3+4*invJ.xelem(6)*Zeta2;
          B.xelem(100) = 4*invJ.xelem(13)*Zeta3+4*invJ.xelem(14)*Zeta2;
          B.xelem(101) = 0;
          B.xelem(102) = 0;
          B.xelem(103) = 0;
          B.xelem(104) = 4*invJ.xelem(13)*Zeta3+4*invJ.xelem(14)*Zeta2;
          B.xelem(105) = 0;
          B.xelem(106) = 4*invJ.xelem(9)*Zeta3+4*invJ.xelem(10)*Zeta2;
          B.xelem(107) = 4*invJ.xelem(5)*Zeta3+4*invJ.xelem(6)*Zeta2;
          B.xelem(108) = 4*invJ.xelem(4)*Zeta3+4*invJ.xelem(6)*Zeta1;
          B.xelem(109) = 0;
          B.xelem(110) = 0;
          B.xelem(111) = 4*invJ.xelem(8)*Zeta3+4*invJ.xelem(10)*Zeta1;
          B.xelem(112) = 0;
          B.xelem(113) = 4*invJ.xelem(12)*Zeta3+4*invJ.xelem(14)*Zeta1;
          B.xelem(114) = 0;
          B.xelem(115) = 4*invJ.xelem(8)*Zeta3+4*invJ.xelem(10)*Zeta1;
          B.xelem(116) = 0;
          B.xelem(117) = 4*invJ.xelem(4)*Zeta3+4*invJ.xelem(6)*Zeta1;
          B.xelem(118) = 4*invJ.xelem(12)*Zeta3+4*invJ.xelem(14)*Zeta1;
          B.xelem(119) = 0;
          B.xelem(120) = 0;
          B.xelem(121) = 0;
          B.xelem(122) = 4*invJ.xelem(12)*Zeta3+4*invJ.xelem(14)*Zeta1;
          B.xelem(123) = 0;
          B.xelem(124) = 4*invJ.xelem(8)*Zeta3+4*invJ.xelem(10)*Zeta1;
          B.xelem(125) = 4*invJ.xelem(4)*Zeta3+4*invJ.xelem(6)*Zeta1;
          B.xelem(126) = 4*invJ.xelem(4)*Zeta4+4*invJ.xelem(7)*Zeta1;
          B.xelem(127) = 0;
          B.xelem(128) = 0;
          B.xelem(129) = 4*invJ.xelem(8)*Zeta4+4*invJ.xelem(11)*Zeta1;
          B.xelem(130) = 0;
          B.xelem(131) = 4*invJ.xelem(12)*Zeta4+4*invJ.xelem(15)*Zeta1;
          B.xelem(132) = 0;
          B.xelem(133) = 4*invJ.xelem(8)*Zeta4+4*invJ.xelem(11)*Zeta1;
          B.xelem(134) = 0;
          B.xelem(135) = 4*invJ.xelem(4)*Zeta4+4*invJ.xelem(7)*Zeta1;
          B.xelem(136) = 4*invJ.xelem(12)*Zeta4+4*invJ.xelem(15)*Zeta1;
          B.xelem(137) = 0;
          B.xelem(138) = 0;
          B.xelem(139) = 0;
          B.xelem(140) = 4*invJ.xelem(12)*Zeta4+4*invJ.xelem(15)*Zeta1;
          B.xelem(141) = 0;
          B.xelem(142) = 4*invJ.xelem(8)*Zeta4+4*invJ.xelem(11)*Zeta1;
          B.xelem(143) = 4*invJ.xelem(4)*Zeta4+4*invJ.xelem(7)*Zeta1;
          B.xelem(144) = 4*invJ.xelem(5)*Zeta4+4*invJ.xelem(7)*Zeta2;
          B.xelem(145) = 0;
          B.xelem(146) = 0;
          B.xelem(147) = 4*invJ.xelem(9)*Zeta4+4*invJ.xelem(11)*Zeta2;
          B.xelem(148) = 0;
          B.xelem(149) = 4*invJ.xelem(13)*Zeta4+4*invJ.xelem(15)*Zeta2;
          B.xelem(150) = 0;
          B.xelem(151) = 4*invJ.xelem(9)*Zeta4+4*invJ.xelem(11)*Zeta2;
          B.xelem(152) = 0;
          B.xelem(153) = 4*invJ.xelem(5)*Zeta4+4*invJ.xelem(7)*Zeta2;
          B.xelem(154) = 4*invJ.xelem(13)*Zeta4+4*invJ.xelem(15)*Zeta2;
          B.xelem(155) = 0;
          B.xelem(156) = 0;
          B.xelem(157) = 0;
          B.xelem(158) = 4*invJ.xelem(13)*Zeta4+4*invJ.xelem(15)*Zeta2;
          B.xelem(159) = 0;
          B.xelem(160) = 4*invJ.xelem(9)*Zeta4+4*invJ.xelem(11)*Zeta2;
          B.xelem(161) = 4*invJ.xelem(5)*Zeta4+4*invJ.xelem(7)*Zeta2;
          B.xelem(162) = 4*invJ.xelem(6)*Zeta4+4*invJ.xelem(7)*Zeta3;
          B.xelem(163) = 0;
          B.xelem(164) = 0;
          B.xelem(165) = 4*invJ.xelem(10)*Zeta4+4*invJ.xelem(11)*Zeta3;
          B.xelem(166) = 0;
          B.xelem(167) = 4*invJ.xelem(14)*Zeta4+4*invJ.xelem(15)*Zeta3;
          B.xelem(168) = 0;
          B.xelem(169) = 4*invJ.xelem(10)*Zeta4+4*invJ.xelem(11)*Zeta3;
          B.xelem(170) = 0;
          B.xelem(171) = 4*invJ.xelem(6)*Zeta4+4*invJ.xelem(7)*Zeta3;
          B.xelem(172) = 4*invJ.xelem(14)*Zeta4+4*invJ.xelem(15)*Zeta3;
          B.xelem(173) = 0;
          B.xelem(174) = 0;
          B.xelem(175) = 0;
          B.xelem(176) = 4*invJ.xelem(14)*Zeta4+4*invJ.xelem(15)*Zeta3;
          B.xelem(177) = 0;
          B.xelem(178) = 4*invJ.xelem(10)*Zeta4+4*invJ.xelem(11)*Zeta3;
          B.xelem(179) = 4*invJ.xelem(6)*Zeta4+4*invJ.xelem(7)*Zeta3;
     }

     virtual Matrix InterpGaussToNodal(FemMatrixType eMatType, const Matrix& taug) const override final {
          return InterpGaussToNodalTpl<double>(eMatType, taug);
     }

     virtual ComplexMatrix InterpGaussToNodal(FemMatrixType eMatType, const ComplexMatrix& taug) const override final {
          return InterpGaussToNodalTpl<std::complex<double> >(eMatType, taug);
     }

     void ScalarInterpMatrix(const ColumnVector& rv, Matrix& Hs, octave_idx_type irow) const override final {
          FEM_ASSERT(rv.numel() == 4);
          FEM_ASSERT(Hs.columns() == 10);
          FEM_ASSERT(irow >= 0);
          FEM_ASSERT(irow < Hs.rows());

          const double Zeta1 = rv.xelem(0);
          const double Zeta2 = rv.xelem(1);
          const double Zeta3 = rv.xelem(2);
          const double Zeta4 = rv.xelem(3);

          const octave_idx_type nrows = Hs.rows();

          Hs.xelem(irow + nrows * 0) = Zeta1*(2*Zeta1-1);
          Hs.xelem(irow + nrows * 1) = Zeta2*(2*Zeta2-1);
          Hs.xelem(irow + nrows * 2) = Zeta3*(2*Zeta3-1);
          Hs.xelem(irow + nrows * 3) = Zeta4*(2*Zeta4-1);
          Hs.xelem(irow + nrows * 4) = 4*Zeta1*Zeta2;
          Hs.xelem(irow + nrows * 5) = 4*Zeta2*Zeta3;
          Hs.xelem(irow + nrows * 6) = 4*Zeta1*Zeta3;
          Hs.xelem(irow + nrows * 7) = 4*Zeta1*Zeta4;
          Hs.xelem(irow + nrows * 8) = 4*Zeta2*Zeta4;
          Hs.xelem(irow + nrows * 9) = 4*Zeta3*Zeta4;
     }

private:
     template <typename T>
     typename PostProcTypeTraits<T>::MatrixType
     InterpGaussToNodalTpl(FemMatrixType eMatType,
                           const typename PostProcTypeTraits<T>::MatrixType& taug) const {
          const IntegrationRule& oIntegRule = GetIntegrationRule(eMatType);
          const octave_idx_type iNumGauss = oIntegRule.iGetNumEvalPoints();
          const octave_idx_type iNumDir = oIntegRule.iGetNumDirections();
          const octave_idx_type iNumNodesCorner = 4;
          const octave_idx_type iNumNodes = nodes.numel();

          ColumnVector rv(iNumDir);

          Matrix H(iNumGauss, iNumNodesCorner);

          for (octave_idx_type i = 0; i < iNumGauss; ++i) {
               for (octave_idx_type j = 0; j < iNumDir; ++j) {
                    rv.xelem(j) = oIntegRule.dGetPosition(i, j);
               }

               ScalarInterpMatrixCornerNodes(rv, H, i);
          }

          typename PostProcTypeTraits<T>::MatrixType taun = H.solve(taug);

          taun.resize(iNumNodes, taun.columns());

          static constexpr struct {
               octave_idx_type ico1, ico2, imid;
          } idxint[] = {
               {0, 1, 4},
               {1, 2, 5},
               {2, 0, 6},
               {0, 3, 7},
               {1, 3, 8},
               {2, 3, 9}
          };

          for (octave_idx_type j = 0; j < taun.columns(); ++j) {
               for (const auto& idx:idxint) {
                    taun.xelem(idx.imid + iNumNodes * j) = 0.5 * (taun.xelem(idx.ico1 + iNumNodes * j) + taun.xelem(idx.ico2 + iNumNodes * j));
               }
          }

          return taun;
     }

     void ScalarInterpMatrixCornerNodes(const ColumnVector& rv, Matrix& Hs, octave_idx_type irow) const {
          FEM_ASSERT(rv.numel() == 4);
          FEM_ASSERT(Hs.columns() == 4);
          FEM_ASSERT(irow >= 0);
          FEM_ASSERT(irow < Hs.rows());

          const double Zeta1 = rv.xelem(0);
          const double Zeta2 = rv.xelem(1);
          const double Zeta3 = rv.xelem(2);
          const double Zeta4 = rv.xelem(3);

          const octave_idx_type nrows = Hs.rows();

          Hs.xelem(irow + nrows * 0) = Zeta1;
          Hs.xelem(irow + nrows * 1) = Zeta2;
          Hs.xelem(irow + nrows * 2) = Zeta3;
          Hs.xelem(irow + nrows * 3) = Zeta4;
     }

     static IntegrationRule& SelectIntegrationRule(FemMatrixType eMatType) {
          switch (eMatType & ~MAT_COLL_PNT_OUTPUT) {
          case MAT_STIFFNESS:
          case MAT_STIFFNESS_SYM:
          case MAT_STIFFNESS_SYM_L:
          case MAT_STIFFNESS_TAU0:
          case VEC_STRESS_CAUCH:
          case VEC_STRAIN_TOTAL:
          case VEC_LOAD_CONSISTENT:
          case VEC_LOAD_LUMPED:
          case VEC_LOAD_FLUID_STRUCT:
          case MAT_THERMAL_COND:
          case MAT_STIFFNESS_ACOUSTICS_RE:
          case MAT_STIFFNESS_ACOUSTICS_IM:
          case MAT_DAMPING_ACOUSTICS_RE:
          case MAT_DAMPING_ACOUSTICS_IM:
          case MAT_STIFFNESS_FLUID_STRUCT_RE:
          case MAT_STIFFNESS_FLUID_STRUCT_IM:
          case MAT_DAMPING_FLUID_STRUCT_RE:
          case MAT_DAMPING_FLUID_STRUCT_IM:
          case VEC_PARTICLE_VELOCITY:
          case VEC_PARTICLE_VELOCITY_C:
               return oIntegStiff;
          case MAT_MASS:
          case MAT_MASS_SYM:
          case MAT_MASS_SYM_L:
          case MAT_STIFFNESS_OMEGA:
          case MAT_STIFFNESS_OMEGA_DOT:
          case MAT_DAMPING_OMEGA:
          case MAT_MASS_FLUID_STRUCT_RE:
          case MAT_MASS_FLUID_STRUCT_IM:
          case VEC_INERTIA_M1:
          case MAT_INERTIA_J:
          case MAT_INERTIA_INV3:
          case MAT_INERTIA_INV4:
          case MAT_INERTIA_INV5:
          case MAT_INERTIA_INV8:
          case MAT_INERTIA_INV9:
          case MAT_HEAT_CAPACITY:
          case MAT_MASS_ACOUSTICS_RE:
          case MAT_MASS_ACOUSTICS_IM:
               return oIntegMass;
          case MAT_MASS_LUMPED:
               return oIntegMassDiag;
          default:
               throw std::runtime_error("tet10: unkown matrix type");
          }
     }

     static IntegrationRule oIntegStiff, oIntegMass, oIntegMassDiag;
};

IntegrationRule Tet10::oIntegStiff;
IntegrationRule Tet10::oIntegMass;
IntegrationRule Tet10::oIntegMassDiag;

class Tet20: public Element3D
{
public:
     Tet20(ElementTypes::TypeId eltype, octave_idx_type id, const Matrix& X, const Material* material, const int32NDArray& nodes, const ElementData& data)
          :Element3D(eltype, id, X, material, nodes, data) {
          FEM_ASSERT(nodes.numel() == 20);
     }

     static void AllocIntegrationRule(FemMatrixType eMatType) {
          switch (eMatType) {
          case MAT_DAMPING:
          case MAT_DAMPING_SYM:
          case MAT_DAMPING_SYM_L:
          case MAT_DAMPING_FLUID_STRUCT_RE:
          case MAT_DAMPING_FLUID_STRUCT_IM:
               AllocIntegrationRule(MAT_MASS);
               AllocIntegrationRule(MAT_STIFFNESS);
               return;
          default:
               break;
          }

          constexpr double a1 = (1 + sqrt(5./14.)) / 4.;
          constexpr double b1 = (1 - sqrt(5./14.)) / 4.;
          constexpr double c1 = 1. / 4.;
          constexpr double d1 = 1. / 14.;
          constexpr double e1 = 11. / 14.;
          constexpr double w1 = -74. / 5625;
          constexpr double w2 = 343. / 45000.;
          constexpr double w3 = 56. / 2250.;

          constexpr octave_idx_type N1 = 11;
          static constexpr double ri1[N1] = {c1, e1, d1, d1, d1, a1, a1, a1, b1, b1, b1};
          static constexpr double si1[N1] = {c1, d1, e1, d1, d1, a1, b1, b1, a1, a1, b1};
          static constexpr double ti1[N1] = {c1, d1, d1, e1, d1, b1, a1, b1, a1, b1, a1};
          static constexpr double wi1[N1] = {w1, w2, w2, w2, w2, w3, w3, w3, w3, w3, w3};

          constexpr octave_idx_type N2 = 27;
          static constexpr double ri2[N2] = {1.431498841332000e-03, 1.127016653792600e-02, 6.350832689629000e-03, 6.350832689629000e-03, 5.000000000000000e-02, 2.817541634481500e-02, 1.127016653792600e-02, 8.872983346207400e-02, 5.000000000000000e-02, 6.350832689629000e-03, 5.000000000000000e-02, 2.817541634481500e-02, 2.817541634481500e-02, 2.218245836551850e-01, 1.250000000000000e-01, 5.000000000000000e-02, 3.936491673103710e-01, 2.218245836551850e-01, 1.127016653792600e-02, 8.872983346207400e-02, 5.000000000000000e-02, 5.000000000000000e-02, 3.936491673103710e-01, 2.218245836551850e-01, 8.872983346207400e-02, 6.985685011586670e-01, 3.936491673103710e-01};
          static constexpr double si2[N2] = {1.127016653792600e-02, 1.431498841332000e-03, 6.350832689629000e-03, 5.000000000000000e-02, 6.350832689629000e-03, 2.817541634481500e-02, 8.872983346207400e-02, 1.127016653792600e-02, 5.000000000000000e-02, 5.000000000000000e-02, 6.350832689629000e-03, 2.817541634481500e-02, 2.218245836551850e-01, 2.817541634481500e-02, 1.250000000000000e-01, 3.936491673103710e-01, 5.000000000000000e-02, 2.218245836551850e-01, 8.872983346207400e-02, 1.127016653792600e-02, 5.000000000000000e-02, 3.936491673103710e-01, 5.000000000000000e-02, 2.218245836551850e-01, 6.985685011586670e-01, 8.872983346207400e-02, 3.936491673103710e-01};
          static constexpr double ti2[N2] = {1.000000000000000e-01, 1.000000000000000e-01, 1.000000000000000e-01, 5.635083268962900e-02, 5.635083268962900e-02, 5.635083268962900e-02, 1.270166537925800e-02, 1.270166537925800e-02, 1.270166537925800e-02, 4.436491673103710e-01, 4.436491673103710e-01, 4.436491673103710e-01, 2.500000000000000e-01, 2.500000000000000e-01, 2.500000000000000e-01, 5.635083268962900e-02, 5.635083268962900e-02, 5.635083268962900e-02, 7.872983346207409e-01, 7.872983346207409e-01, 7.872983346207409e-01, 4.436491673103710e-01, 4.436491673103710e-01, 4.436491673103710e-01, 1.000000000000000e-01, 1.000000000000000e-01, 1.000000000000000e-01};
          static constexpr double wi2[N2] = {3.068198819728420e-05, 3.068198819728420e-05, 4.909118111565470e-05, 2.177926162424280e-04, 2.177926162424280e-04, 3.484681859878840e-04, 2.415587821057510e-04, 2.415587821057510e-04, 3.864940513692010e-04, 9.662351284230000e-04, 9.662351284230000e-04, 1.545976205477000e-03, 6.858710562414000e-03, 6.858710562414000e-03, 1.097393689986300e-02, 7.607153074595000e-03, 7.607153074595000e-03, 1.217144491935200e-02, 1.901788268649000e-03, 1.901788268649000e-03, 3.042861229838000e-03, 1.349962850858600e-02, 1.349962850858600e-02, 2.159940561373800e-02, 1.497274736708400e-02, 1.497274736708400e-02, 2.395639578733400e-02};

          constexpr octave_idx_type N3 = 64;
          static constexpr double ri3[N3] = {3.347157145940000e-04, 4.486065274832000e-03, 1.590903418873000e-03, 3.229877570553000e-03, 1.590903418873000e-03, 2.132226325753900e-02, 7.561562178966000e-03, 1.535160449744700e-02, 3.229877570553000e-03, 4.328879995600800e-02, 1.535160449744700e-02, 3.116707302911400e-02, 4.486065274832000e-03, 6.012499793871600e-02, 2.132226325753900e-02, 4.328879995600800e-02, 1.590903418873000e-03, 2.132226325753900e-02, 7.561562178966000e-03, 1.535160449744700e-02, 7.561562178966000e-03, 1.013446935278680e-01, 3.594009661935300e-02, 7.296615908748100e-02, 1.535160449744700e-02, 2.057516180032910e-01, 7.296615908748100e-02, 1.481370634132570e-01, 2.132226325753900e-02, 2.857740482736200e-01, 1.013446935278680e-01, 2.057516180032910e-01, 3.229877570553000e-03, 4.328879995600800e-02, 1.535160449744700e-02, 3.116707302911400e-02, 1.535160449744700e-02, 2.057516180032910e-01, 7.296615908748100e-02, 1.481370634132570e-01, 3.116707302911400e-02, 4.177202262625760e-01, 1.481370634132570e-01, 3.007502358784330e-01, 4.328879995600800e-02, 5.801830443098590e-01, 2.057516180032910e-01, 4.177202262625760e-01, 4.486065274832000e-03, 6.012499793871600e-02, 2.132226325753900e-02, 4.328879995600800e-02, 2.132226325753900e-02, 2.857740482736200e-01, 1.013446935278680e-01, 2.057516180032910e-01, 4.328879995600800e-02, 5.801830443098590e-01, 2.057516180032910e-01, 4.177202262625760e-01, 6.012499793871600e-02, 8.058320946447630e-01, 2.857740482736200e-01, 5.801830443098590e-01};
          static constexpr double si3[N3] = {4.486065274832000e-03, 3.347157145940000e-04, 3.229877570553000e-03, 1.590903418873000e-03, 2.132226325753900e-02, 1.590903418873000e-03, 1.535160449744700e-02, 7.561562178966000e-03, 4.328879995600800e-02, 3.229877570553000e-03, 3.116707302911400e-02, 1.535160449744700e-02, 6.012499793871600e-02, 4.486065274832000e-03, 4.328879995600800e-02, 2.132226325753900e-02, 2.132226325753900e-02, 1.590903418873000e-03, 1.535160449744700e-02, 7.561562178966000e-03, 1.013446935278680e-01, 7.561562178966000e-03, 7.296615908748100e-02, 3.594009661935300e-02, 2.057516180032910e-01, 1.535160449744700e-02, 1.481370634132570e-01, 7.296615908748100e-02, 2.857740482736200e-01, 2.132226325753900e-02, 2.057516180032910e-01, 1.013446935278680e-01, 4.328879995600800e-02, 3.229877570553000e-03, 3.116707302911400e-02, 1.535160449744700e-02, 2.057516180032910e-01, 1.535160449744700e-02, 1.481370634132570e-01, 7.296615908748100e-02, 4.177202262625760e-01, 3.116707302911400e-02, 3.007502358784330e-01, 1.481370634132570e-01, 5.801830443098590e-01, 4.328879995600800e-02, 4.177202262625760e-01, 2.057516180032910e-01, 6.012499793871600e-02, 4.486065274832000e-03, 4.328879995600800e-02, 2.132226325753900e-02, 2.857740482736200e-01, 2.132226325753900e-02, 2.057516180032910e-01, 1.013446935278680e-01, 5.801830443098590e-01, 4.328879995600800e-02, 4.177202262625760e-01, 2.057516180032910e-01, 8.058320946447630e-01, 6.012499793871600e-02, 5.801830443098590e-01, 2.857740482736200e-01};
          static constexpr double ti3[N3] = {6.461106321354800e-02, 6.461106321354800e-02, 6.461106321354800e-02, 6.461106321354800e-02, 4.651867752656100e-02, 4.651867752656100e-02, 4.651867752656100e-02, 4.651867752656100e-02, 2.291316667641300e-02, 2.291316667641300e-02, 2.291316667641300e-02, 2.291316667641300e-02, 4.820780989426000e-03, 4.820780989426000e-03, 4.820780989426000e-03, 4.820780989426000e-03, 3.070963115311590e-01, 3.070963115311590e-01, 3.070963115311590e-01, 3.070963115311590e-01, 2.211032225007380e-01, 2.211032225007380e-01, 2.211032225007380e-01, 2.211032225007380e-01, 1.089062557068340e-01, 1.089062557068340e-01, 1.089062557068340e-01, 1.089062557068340e-01, 2.291316667641300e-02, 2.291316667641300e-02, 2.291316667641300e-02, 2.291316667641300e-02, 6.234718442658670e-01, 6.234718442658670e-01, 6.234718442658670e-01, 6.234718442658670e-01, 4.488872992916900e-01, 4.488872992916900e-01, 4.488872992916900e-01, 4.488872992916900e-01, 2.211032225007380e-01, 2.211032225007380e-01, 2.211032225007380e-01, 2.211032225007380e-01, 4.651867752656100e-02, 4.651867752656100e-02, 4.651867752656100e-02, 4.651867752656100e-02, 8.659570925834790e-01, 8.659570925834790e-01, 8.659570925834790e-01, 8.659570925834790e-01, 6.234718442658670e-01, 6.234718442658670e-01, 6.234718442658670e-01, 6.234718442658670e-01, 3.070963115311590e-01, 3.070963115311590e-01, 3.070963115311590e-01, 3.070963115311590e-01, 6.461106321354800e-02, 6.461106321354800e-02, 6.461106321354800e-02, 6.461106321354800e-02};
          static constexpr double wi3[N3] = {1.761084870822600e-06, 1.761084870822600e-06, 3.301615549885100e-06, 3.301615549885100e-06, 1.569257503335800e-05, 1.569257503335800e-05, 2.941984830275260e-05, 2.941984830275260e-05, 3.185931686560010e-05, 3.185931686560010e-05, 5.972864665122530e-05, 5.972864665122530e-05, 2.360313944207820e-05, 2.360313944207820e-05, 4.425027634907300e-05, 4.425027634907300e-05, 7.458679166511140e-05, 7.458679166511140e-05, 1.398325062338110e-04, 1.398325062338110e-04, 6.646237464725280e-04, 6.646237464725280e-04, 1.246011553748000e-03, 1.246011553748000e-03, 1.349329762022000e-03, 1.349329762022000e-03, 2.529672588768000e-03, 2.529672588768000e-03, 9.996579230088670e-04, 9.996579230088670e-04, 1.874121002261000e-03, 1.874121002261000e-03, 3.074301219528830e-04, 3.074301219528830e-04, 5.763584072291740e-04, 5.763584072291740e-04, 2.739430867978000e-03, 2.739430867978000e-03, 5.135781756688000e-03, 5.135781756688000e-03, 5.561636370626000e-03, 5.561636370626000e-03, 1.042674627912200e-02, 1.042674627912200e-02, 4.120367029080000e-03, 4.120367029080000e-03, 7.724708831375000e-03, 7.724708831375000e-03, 3.163437496694650e-04, 3.163437496694650e-04, 5.930693405649470e-04, 5.930693405649470e-04, 2.818857915521000e-03, 2.818857915521000e-03, 5.284688592238000e-03, 5.284688592238000e-03, 5.722890433136000e-03, 5.722890433136000e-03, 1.072905931870600e-02, 1.072905931870600e-02, 4.239832934111000e-03, 4.239832934111000e-03, 7.948679008091999e-03, 7.948679008091999e-03};

          constexpr double a4 = (5. + 3. * sqrt(5)) / 20.;
          constexpr double b4 = (5. - sqrt(5.)) / 20.;
          constexpr double c4 = 1. / 24.;
          constexpr octave_idx_type N4 = 4;

          static constexpr double ri4[N4] = {b4, a4, b4, b4};
          static constexpr double si4[N4] = {b4, b4, a4, b4};
          static constexpr double ti4[N4] = {b4, b4, b4, a4};
          static constexpr double wi4[N4] = {c4, c4, c4, c4};

          constexpr octave_idx_type N5 = 20;

          static constexpr double ri5[N5] = {0.0000000000000000e+00, 3.3333333333333331e-01, 6.6666666666666663e-01, 1.0000000000000000e+00, 6.6666666666666663e-01, 3.3333333333333331e-01, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 3.3333333333333331e-01, 6.6666666666666663e-01, 3.3333333333333331e-01, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 3.3333333333333331e-01, 3.3333333333333331e-01, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00};
          static constexpr double si5[N5] = {0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 3.3333333333333331e-01, 6.6666666666666663e-01, 1.0000000000000000e+00, 6.6666666666666663e-01, 3.3333333333333331e-01, 3.3333333333333331e-01, 0.0000000000000000e+00, 3.3333333333333331e-01, 6.6666666666666663e-01, 3.3333333333333331e-01, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 3.3333333333333331e-01, 0.0000000000000000e+00, 0.0000000000000000e+00};
          static constexpr double ti5[N5] = {0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 3.3333333333333331e-01, 3.3333333333333331e-01, 3.3333333333333331e-01, 3.3333333333333331e-01, 3.3333333333333331e-01, 3.3333333333333331e-01, 6.6666666666666663e-01, 6.6666666666666663e-01, 6.6666666666666663e-01, 1.0000000000000000e+00};
          static constexpr double wi5[N5] = {8.333333333333333e-03, 8.333333333333333e-03, 8.333333333333333e-03, 8.333333333333333e-03, 8.333333333333333e-03, 8.333333333333333e-03, 8.333333333333333e-03, 8.333333333333333e-03, 8.333333333333333e-03, 8.333333333333333e-03, 8.333333333333333e-03, 8.333333333333333e-03, 8.333333333333333e-03, 8.333333333333333e-03, 8.333333333333333e-03, 8.333333333333333e-03, 8.333333333333333e-03, 8.333333333333333e-03, 8.333333333333333e-03, 8.333333333333333e-03};

          static constexpr octave_idx_type M = 5;
          static constexpr octave_idx_type N[M] = {N1, N2, N3, N4, N5};
          static constexpr const double* r[M] = {ri1, ri2, ri3, ri4, ri5};
          static constexpr const double* s[M] = {si1, si2, si3, si4, si5};
          static constexpr const double* t[M] = {ti1, ti2, ti3, ti4, ti5};
          static constexpr const double* w[M] = {wi1, wi2, wi3, wi4, wi5};

          const octave_idx_type iIntegRule = SelectIntegrationRule(eMatType);

          FEM_ASSERT(iIntegRule >= 0);
          FEM_ASSERT(iIntegRule < M);

          if (!rgIntegRule[iIntegRule].iGetNumEvalPoints()) {
               rgIntegRule[iIntegRule].SetNumEvalPoints(N[iIntegRule], 3);

               for (octave_idx_type i = 0; i < N[iIntegRule]; ++i) {
                    rgIntegRule[iIntegRule].SetPosition(i, 0, r[iIntegRule][i]);
                    rgIntegRule[iIntegRule].SetPosition(i, 1, s[iIntegRule][i]);
                    rgIntegRule[iIntegRule].SetPosition(i, 2, t[iIntegRule][i]);
                    rgIntegRule[iIntegRule].SetWeight(i, w[iIntegRule][i]);
               }
          }
     }

     virtual const IntegrationRule& GetIntegrationRule(FemMatrixType eMatType) const override final {
          const octave_idx_type iIntegRule = SelectIntegrationRule(eMatType);

          FEM_ASSERT(iIntegRule >= 0);
          FEM_ASSERT(iIntegRule < RNUM);
          FEM_ASSERT(rgIntegRule[iIntegRule].iGetNumEvalPoints() > 0);

          return rgIntegRule[iIntegRule];
     }

protected:
     virtual Matrix InterpGaussToNodal(FemMatrixType eMatType, const Matrix& taug) const override final {
          return InterpGaussToNodalTpl<double>(eMatType, taug);
     }

     virtual ComplexMatrix InterpGaussToNodal(FemMatrixType eMatType, const ComplexMatrix& taug) const override final {
          return InterpGaussToNodalTpl<std::complex<double> >(eMatType, taug);
     }

     void ScalarInterpMatrix(const ColumnVector& rv, Matrix& Hs, octave_idx_type irow) const override final {
          FEM_ASSERT(rv.numel() == 3);
          FEM_ASSERT(Hs.columns() == 20);
          FEM_ASSERT(irow >= 0);
          FEM_ASSERT(irow < Hs.rows());

          const double r = rv.xelem(0);
          const double s = rv.xelem(1);
          const double t = rv.xelem(2);
          const octave_idx_type nrows = Hs.rows();

          Hs.xelem(irow) = ((3*((-t)-s-r+1)-2)*(3*((-t)-s-r+1)-1)*((-t)-s-r+1))/2.0E+0;
          Hs.xelem(nrows+irow) = (9.0E+0*r*(3*((-t)-s-r+1)-1)*((-t)-s-r+1))/2.0E+0;
          Hs.xelem(2*nrows+irow) = (9.0E+0*r*(3*r-1)*((-t)-s-r+1))/2.0E+0;
          Hs.xelem(3*nrows+irow) = (r*(3*r-2)*(3*r-1))/2.0E+0;
          Hs.xelem(4*nrows+irow) = (9.0E+0*r*(3*r-1)*s)/2.0E+0;
          Hs.xelem(5*nrows+irow) = (9.0E+0*r*s*(3*s-1))/2.0E+0;
          Hs.xelem(6*nrows+irow) = (s*(3*s-2)*(3*s-1))/2.0E+0;
          Hs.xelem(7*nrows+irow) = (9.0E+0*s*(3*s-1)*((-t)-s-r+1))/2.0E+0;
          Hs.xelem(8*nrows+irow) = (9.0E+0*s*(3*((-t)-s-r+1)-1)*((-t)-s-r+1))/2.0E+0;
          Hs.xelem(9*nrows+irow) = 27*r*s*((-t)-s-r+1);
          Hs.xelem(10*nrows+irow) = (9.0E+0*r*(3*r-1)*t)/2.0E+0;
          Hs.xelem(11*nrows+irow) = 27*r*s*t;
          Hs.xelem(12*nrows+irow) = (9.0E+0*s*(3*s-1)*t)/2.0E+0;
          Hs.xelem(13*nrows+irow) = 27*s*((-t)-s-r+1)*t;
          Hs.xelem(14*nrows+irow) = (9.0E+0*(3*((-t)-s-r+1)-1)*((-t)-s-r+1)*t)/2.0E+0;
          Hs.xelem(15*nrows+irow) = 27*r*((-t)-s-r+1)*t;
          Hs.xelem(16*nrows+irow) = (9.0E+0*r*t*(3*t-1))/2.0E+0;
          Hs.xelem(17*nrows+irow) = (9.0E+0*s*t*(3*t-1))/2.0E+0;
          Hs.xelem(18*nrows+irow) = (9.0E+0*((-t)-s-r+1)*t*(3*t-1))/2.0E+0;
          Hs.xelem(19*nrows+irow) = (t*(3*t-2)*(3*t-1))/2.0E+0;
     }

     virtual void ScalarInterpMatrixDer(const ColumnVector& rv, Matrix& Hd) const override final {
          const double r = rv.xelem(0);
          const double s = rv.xelem(1);
          const double t = rv.xelem(2);

          FEM_ASSERT(rv.numel() == 3);
          FEM_ASSERT(Hd.rows() == 20);
          FEM_ASSERT(Hd.columns() == 3);

          Hd.xelem(0) = ((-3.0E+0)*(3*((-t)-s-r+1)-1)*((-t)-s-r+1))/2.0E+0+((-3.0E+0)*(3*((-t)-s-r+1)-2)*((-t)-s-r+1))/2.0E+0-((3*((-t)-s-r+1)-2)*(3*((-t)-s-r+1)-1))/2.0E+0;
          Hd.xelem(1) = (9.0E+0*(3*((-t)-s-r+1)-1)*((-t)-s-r+1))/2.0E+0+((-2.7E+1)*r*((-t)-s-r+1))/2.0E+0+((-9.0E+0)*r*(3*((-t)-s-r+1)-1))/2.0E+0;
          Hd.xelem(2) = (9.0E+0*(3*r-1)*((-t)-s-r+1))/2.0E+0+(2.7E+1*r*((-t)-s-r+1))/2.0E+0+((-9.0E+0)*r*(3*r-1))/2.0E+0;
          Hd.xelem(3) = ((3*r-2)*(3*r-1))/2.0E+0+(3.0E+0*r*(3*r-1))/2.0E+0+(3.0E+0*r*(3*r-2))/2.0E+0;
          Hd.xelem(4) = (9.0E+0*(3*r-1)*s)/2.0E+0+(2.7E+1*r*s)/2.0E+0;
          Hd.xelem(5) = (9.0E+0*s*(3*s-1))/2.0E+0;
          Hd.xelem(6) = 0;
          Hd.xelem(7) = ((-9.0E+0)*s*(3*s-1))/2.0E+0;
          Hd.xelem(8) = ((-2.7E+1)*s*((-t)-s-r+1))/2.0E+0+((-9.0E+0)*s*(3*((-t)-s-r+1)-1))/2.0E+0;
          Hd.xelem(9) = 27*s*((-t)-s-r+1)-27*r*s;
          Hd.xelem(10) = (9.0E+0*(3*r-1)*t)/2.0E+0+(2.7E+1*r*t)/2.0E+0;
          Hd.xelem(11) = 27*s*t;
          Hd.xelem(12) = 0;
          Hd.xelem(13) = -27*s*t;
          Hd.xelem(14) = ((-2.7E+1)*((-t)-s-r+1)*t)/2.0E+0+((-9.0E+0)*(3*((-t)-s-r+1)-1)*t)/2.0E+0;
          Hd.xelem(15) = 27*((-t)-s-r+1)*t-27*r*t;
          Hd.xelem(16) = (9.0E+0*t*(3*t-1))/2.0E+0;
          Hd.xelem(17) = 0;
          Hd.xelem(18) = ((-9.0E+0)*t*(3*t-1))/2.0E+0;
          Hd.xelem(19) = 0;
          Hd.xelem(20) = ((-3.0E+0)*(3*((-t)-s-r+1)-1)*((-t)-s-r+1))/2.0E+0+((-3.0E+0)*(3*((-t)-s-r+1)-2)*((-t)-s-r+1))/2.0E+0-((3*((-t)-s-r+1)-2)*(3*((-t)-s-r+1)-1))/2.0E+0;
          Hd.xelem(21) = ((-2.7E+1)*r*((-t)-s-r+1))/2.0E+0+((-9.0E+0)*r*(3*((-t)-s-r+1)-1))/2.0E+0;
          Hd.xelem(22) = ((-9.0E+0)*r*(3*r-1))/2.0E+0;
          Hd.xelem(23) = 0;
          Hd.xelem(24) = (9.0E+0*r*(3*r-1))/2.0E+0;
          Hd.xelem(25) = (9.0E+0*r*(3*s-1))/2.0E+0+(2.7E+1*r*s)/2.0E+0;
          Hd.xelem(26) = ((3*s-2)*(3*s-1))/2.0E+0+(3.0E+0*s*(3*s-1))/2.0E+0+(3.0E+0*s*(3*s-2))/2.0E+0;
          Hd.xelem(27) = (9.0E+0*(3*s-1)*((-t)-s-r+1))/2.0E+0+(2.7E+1*s*((-t)-s-r+1))/2.0E+0+((-9.0E+0)*s*(3*s-1))/2.0E+0;
          Hd.xelem(28) = (9.0E+0*(3*((-t)-s-r+1)-1)*((-t)-s-r+1))/2.0E+0+((-2.7E+1)*s*((-t)-s-r+1))/2.0E+0+((-9.0E+0)*s*(3*((-t)-s-r+1)-1))/2.0E+0;
          Hd.xelem(29) = 27*r*((-t)-s-r+1)-27*r*s;
          Hd.xelem(30) = 0;
          Hd.xelem(31) = 27*r*t;
          Hd.xelem(32) = (9.0E+0*(3*s-1)*t)/2.0E+0+(2.7E+1*s*t)/2.0E+0;
          Hd.xelem(33) = 27*((-t)-s-r+1)*t-27*s*t;
          Hd.xelem(34) = ((-2.7E+1)*((-t)-s-r+1)*t)/2.0E+0+((-9.0E+0)*(3*((-t)-s-r+1)-1)*t)/2.0E+0;
          Hd.xelem(35) = -27*r*t;
          Hd.xelem(36) = 0;
          Hd.xelem(37) = (9.0E+0*t*(3*t-1))/2.0E+0;
          Hd.xelem(38) = ((-9.0E+0)*t*(3*t-1))/2.0E+0;
          Hd.xelem(39) = 0;
          Hd.xelem(40) = ((-3.0E+0)*(3*((-t)-s-r+1)-1)*((-t)-s-r+1))/2.0E+0+((-3.0E+0)*(3*((-t)-s-r+1)-2)*((-t)-s-r+1))/2.0E+0-((3*((-t)-s-r+1)-2)*(3*((-t)-s-r+1)-1))/2.0E+0;
          Hd.xelem(41) = ((-2.7E+1)*r*((-t)-s-r+1))/2.0E+0+((-9.0E+0)*r*(3*((-t)-s-r+1)-1))/2.0E+0;
          Hd.xelem(42) = ((-9.0E+0)*r*(3*r-1))/2.0E+0;
          Hd.xelem(43) = 0;
          Hd.xelem(44) = 0;
          Hd.xelem(45) = 0;
          Hd.xelem(46) = 0;
          Hd.xelem(47) = ((-9.0E+0)*s*(3*s-1))/2.0E+0;
          Hd.xelem(48) = ((-2.7E+1)*s*((-t)-s-r+1))/2.0E+0+((-9.0E+0)*s*(3*((-t)-s-r+1)-1))/2.0E+0;
          Hd.xelem(49) = -27*r*s;
          Hd.xelem(50) = (9.0E+0*r*(3*r-1))/2.0E+0;
          Hd.xelem(51) = 27*r*s;
          Hd.xelem(52) = (9.0E+0*s*(3*s-1))/2.0E+0;
          Hd.xelem(53) = 27*s*((-t)-s-r+1)-27*s*t;
          Hd.xelem(54) = ((-2.7E+1)*((-t)-s-r+1)*t)/2.0E+0+((-9.0E+0)*(3*((-t)-s-r+1)-1)*t)/2.0E+0+(9.0E+0*(3*((-t)-s-r+1)-1)*((-t)-s-r+1))/2.0E+0;
          Hd.xelem(55) = 27*r*((-t)-s-r+1)-27*r*t;
          Hd.xelem(56) = (9.0E+0*r*(3*t-1))/2.0E+0+(2.7E+1*r*t)/2.0E+0;
          Hd.xelem(57) = (9.0E+0*s*(3*t-1))/2.0E+0+(2.7E+1*s*t)/2.0E+0;
          Hd.xelem(58) = ((-9.0E+0)*t*(3*t-1))/2.0E+0+(9.0E+0*((-t)-s-r+1)*(3*t-1))/2.0E+0+(2.7E+1*((-t)-s-r+1)*t)/2.0E+0;
          Hd.xelem(59) = ((3*t-2)*(3*t-1))/2.0E+0+(3.0E+0*t*(3*t-1))/2.0E+0+(3.0E+0*t*(3*t-2))/2.0E+0;
     }
private:
     enum IntegRuleType {
          R1 = 0,
          R2,
          R3,
          R4,
          R5,
          RNUM
     };

     template <typename T>
     typename PostProcTypeTraits<T>::MatrixType
     InterpGaussToNodalTpl(FemMatrixType eMatType,
                           const typename PostProcTypeTraits<T>::MatrixType& taug) const {
          const IntegrationRule& oIntegRule = GetIntegrationRule(eMatType);
          const octave_idx_type iNumGauss = oIntegRule.iGetNumEvalPoints();
          const octave_idx_type iNumDir = oIntegRule.iGetNumDirections();
          const octave_idx_type iNumNodes = nodes.numel();

          ColumnVector rv(iNumDir);

          FEM_ASSERT(iNumGauss >= iNumNodes);

          Matrix H(iNumGauss, iNumNodes);

          for (octave_idx_type i = 0; i < iNumGauss; ++i) {
               for (octave_idx_type j = 0; j < iNumDir; ++j) {
                    rv.xelem(j) = oIntegRule.dGetPosition(i, j);
               }

               ScalarInterpMatrix(rv, H, i);
          }

          octave_idx_type info = -1;

          const auto taun = H.solve(taug, info);

          FEM_ASSERT(info == 0);

          return taun;
     }

     static octave_idx_type SelectIntegrationRule(FemMatrixType eMatType) {
          switch (eMatType & ~MAT_COLL_PNT_OUTPUT) {
          case MAT_STIFFNESS:
          case MAT_STIFFNESS_SYM:
          case MAT_STIFFNESS_SYM_L:
          case MAT_STIFFNESS_TAU0:
          case MAT_STIFFNESS_FLUID_STRUCT_RE:
          case MAT_STIFFNESS_FLUID_STRUCT_IM:
          case MAT_DAMPING_FLUID_STRUCT_RE:
          case MAT_DAMPING_FLUID_STRUCT_IM:
          case VEC_LOAD_CONSISTENT:
          case VEC_LOAD_LUMPED:
          case VEC_LOAD_THERMAL:
          case VEC_LOAD_ACOUSTICS:
          case VEC_LOAD_FLUID_STRUCT:
          case MAT_THERMAL_COND:
          case MAT_STIFFNESS_ACOUSTICS_RE:
          case MAT_STIFFNESS_ACOUSTICS_IM:
          case MAT_DAMPING_ACOUSTICS_RE:
          case MAT_DAMPING_ACOUSTICS_IM:
               return R1;

          case MAT_MASS:
          case MAT_MASS_SYM:
          case MAT_MASS_SYM_L:
          case MAT_STIFFNESS_OMEGA:
          case MAT_STIFFNESS_OMEGA_DOT:
          case MAT_DAMPING_OMEGA:
          case MAT_MASS_FLUID_STRUCT_RE:
          case MAT_MASS_FLUID_STRUCT_IM:
          case VEC_INERTIA_M1:
          case MAT_INERTIA_J:
          case MAT_INERTIA_INV3:
          case MAT_INERTIA_INV4:
          case MAT_INERTIA_INV5:
          case MAT_INERTIA_INV8:
          case MAT_INERTIA_INV9:
          case SCA_TOT_MASS:
          case MAT_HEAT_CAPACITY:
          case MAT_MASS_ACOUSTICS_RE:
          case MAT_MASS_ACOUSTICS_IM:
               return R3;

          case VEC_STRESS_CAUCH:
          case VEC_STRAIN_TOTAL:
          case SCA_STRESS_VMIS:
          case VEC_PARTICLE_VELOCITY:
          case VEC_PARTICLE_VELOCITY_C:
               return R3;

          case MAT_MASS_LUMPED:
               return R5;

          default:
               throw std::runtime_error("tet20: matrix type not supported");
          }
     }

     void ScalarInterpMatrixCornerNodes(const ColumnVector& rv, Matrix& Hs, octave_idx_type irow) const {
          FEM_ASSERT(rv.numel() == 3);
          FEM_ASSERT(Hs.columns() == 4);
          FEM_ASSERT(irow >= 0);
          FEM_ASSERT(irow < Hs.rows());

          const double r = rv.xelem(0);
          const double s = rv.xelem(1);
          const double t = rv.xelem(2);

          const octave_idx_type nrows = Hs.rows();

          Hs.xelem(0*nrows+irow) = 1. - r - s - t;
          Hs.xelem(1*nrows+irow) = r;
          Hs.xelem(2*nrows+irow) = s;
          Hs.xelem(3*nrows+irow) = t;
     }

     static array<IntegrationRule, RNUM> rgIntegRule;
};

array<IntegrationRule, Tet20::RNUM> Tet20::rgIntegRule;

class ShapeTria6 {
public:
     static constexpr octave_idx_type iGetNumNodes() {
          return 6;
     }

     static constexpr octave_idx_type iGetNumDirections() {
          return 3;
     }

     static constexpr octave_idx_type iGetNumDofNode() {
          return 3;
     }

     static constexpr octave_idx_type iGetNumEqualityConstr() {
          return 1;
     }

     static constexpr octave_idx_type iGetNumInequalityConstr() {
          return 0;
     }

     static double EqualityConstr(const ColumnVector& rv) {
          return rv.xelem(0) + rv.xelem(1) + rv.xelem(2) - 1;
     }

     static double InequalityConstr(const ColumnVector& rv) {
          return 0.;
     }

     static void GetElemLimits(ColumnVector& rmin, ColumnVector& rmax) {
          FEM_ASSERT(rmin.rows() == iGetNumDirections());
          FEM_ASSERT(rmax.rows() == rmin.rows());

          for (octave_idx_type i = 0; i < rmin.rows(); ++i) {
               rmin.xelem(i) = 0.;
          }

          for (octave_idx_type i = 0; i < rmax.rows(); ++i) {
               rmax.xelem(i) = 1.;
          }
     }

     static void ScalarInterpMatrix(const ColumnVector& rv, Matrix& HA, octave_idx_type irow) {
          FEM_ASSERT(rv.rows() == 3);
          FEM_ASSERT(HA.rows() > irow);
          FEM_ASSERT(HA.columns() == 6);

          const double Zeta1 = rv.xelem(0);
          const double Zeta2 = rv.xelem(1);
          const double Zeta3 = rv.xelem(2);

          const octave_idx_type nrows = HA.rows();

          HA.xelem(irow + nrows * 0) = Zeta1*(2*Zeta1-1);
          HA.xelem(irow + nrows * 1) = Zeta2*(2*Zeta2-1);
          HA.xelem(irow + nrows * 2) = Zeta3*(2*Zeta3-1);
          HA.xelem(irow + nrows * 3) = 4*Zeta1*Zeta2;
          HA.xelem(irow + nrows * 4) = 4*Zeta2*Zeta3;
          HA.xelem(irow + nrows * 5) = 4*Zeta1*Zeta3;
     }

     static void VectorInterpMatrix(const ColumnVector& rv, Matrix& Hf) {
          FEM_ASSERT(rv.rows() == 3);
          FEM_ASSERT(Hf.rows() == 3);
          FEM_ASSERT(Hf.columns() == 18);

          const double Zeta1 = rv.xelem(0);
          const double Zeta2 = rv.xelem(1);
          const double Zeta3 = rv.xelem(2);

          Hf.xelem(0) = Zeta1*(2*Zeta1-1);
          Hf.xelem(1) = 0;
          Hf.xelem(2) = 0;
          Hf.xelem(3) = 0;
          Hf.xelem(4) = Zeta1*(2*Zeta1-1);
          Hf.xelem(5) = 0;
          Hf.xelem(6) = 0;
          Hf.xelem(7) = 0;
          Hf.xelem(8) = Zeta1*(2*Zeta1-1);
          Hf.xelem(9) = Zeta2*(2*Zeta2-1);
          Hf.xelem(10) = 0;
          Hf.xelem(11) = 0;
          Hf.xelem(12) = 0;
          Hf.xelem(13) = Zeta2*(2*Zeta2-1);
          Hf.xelem(14) = 0;
          Hf.xelem(15) = 0;
          Hf.xelem(16) = 0;
          Hf.xelem(17) = Zeta2*(2*Zeta2-1);
          Hf.xelem(18) = Zeta3*(2*Zeta3-1);
          Hf.xelem(19) = 0;
          Hf.xelem(20) = 0;
          Hf.xelem(21) = 0;
          Hf.xelem(22) = Zeta3*(2*Zeta3-1);
          Hf.xelem(23) = 0;
          Hf.xelem(24) = 0;
          Hf.xelem(25) = 0;
          Hf.xelem(26) = Zeta3*(2*Zeta3-1);
          Hf.xelem(27) = 4*Zeta1*Zeta2;
          Hf.xelem(28) = 0;
          Hf.xelem(29) = 0;
          Hf.xelem(30) = 0;
          Hf.xelem(31) = 4*Zeta1*Zeta2;
          Hf.xelem(32) = 0;
          Hf.xelem(33) = 0;
          Hf.xelem(34) = 0;
          Hf.xelem(35) = 4*Zeta1*Zeta2;
          Hf.xelem(36) = 4*Zeta2*Zeta3;
          Hf.xelem(37) = 0;
          Hf.xelem(38) = 0;
          Hf.xelem(39) = 0;
          Hf.xelem(40) = 4*Zeta2*Zeta3;
          Hf.xelem(41) = 0;
          Hf.xelem(42) = 0;
          Hf.xelem(43) = 0;
          Hf.xelem(44) = 4*Zeta2*Zeta3;
          Hf.xelem(45) = 4*Zeta1*Zeta3;
          Hf.xelem(46) = 0;
          Hf.xelem(47) = 0;
          Hf.xelem(48) = 0;
          Hf.xelem(49) = 4*Zeta1*Zeta3;
          Hf.xelem(50) = 0;
          Hf.xelem(51) = 0;
          Hf.xelem(52) = 0;
          Hf.xelem(53) = 4*Zeta1*Zeta3;
     }

     static void VectorInterpMatrixDerR(const ColumnVector& rv, Matrix& dHf_dr) {
          FEM_ASSERT(rv.rows() == 3);
          FEM_ASSERT(dHf_dr.rows() == 3);
          FEM_ASSERT(dHf_dr.columns() == 18);

          const double Zeta1 = rv.xelem(0);
          const double Zeta2 = rv.xelem(1);
          const double Zeta3 = rv.xelem(2);

          dHf_dr.xelem(0) = 4*Zeta1-1;
          dHf_dr.xelem(1) = 0;
          dHf_dr.xelem(2) = 0;
          dHf_dr.xelem(3) = 0;
          dHf_dr.xelem(4) = 4*Zeta1-1;
          dHf_dr.xelem(5) = 0;
          dHf_dr.xelem(6) = 0;
          dHf_dr.xelem(7) = 0;
          dHf_dr.xelem(8) = 4*Zeta1-1;
          dHf_dr.xelem(9) = 0;
          dHf_dr.xelem(10) = 0;
          dHf_dr.xelem(11) = 0;
          dHf_dr.xelem(12) = 0;
          dHf_dr.xelem(13) = 0;
          dHf_dr.xelem(14) = 0;
          dHf_dr.xelem(15) = 0;
          dHf_dr.xelem(16) = 0;
          dHf_dr.xelem(17) = 0;
          dHf_dr.xelem(18) = 1-4*Zeta3;
          dHf_dr.xelem(19) = 0;
          dHf_dr.xelem(20) = 0;
          dHf_dr.xelem(21) = 0;
          dHf_dr.xelem(22) = 1-4*Zeta3;
          dHf_dr.xelem(23) = 0;
          dHf_dr.xelem(24) = 0;
          dHf_dr.xelem(25) = 0;
          dHf_dr.xelem(26) = 1-4*Zeta3;
          dHf_dr.xelem(27) = 4*Zeta2;
          dHf_dr.xelem(28) = 0;
          dHf_dr.xelem(29) = 0;
          dHf_dr.xelem(30) = 0;
          dHf_dr.xelem(31) = 4*Zeta2;
          dHf_dr.xelem(32) = 0;
          dHf_dr.xelem(33) = 0;
          dHf_dr.xelem(34) = 0;
          dHf_dr.xelem(35) = 4*Zeta2;
          dHf_dr.xelem(36) = -4*Zeta2;
          dHf_dr.xelem(37) = 0;
          dHf_dr.xelem(38) = 0;
          dHf_dr.xelem(39) = 0;
          dHf_dr.xelem(40) = -4*Zeta2;
          dHf_dr.xelem(41) = 0;
          dHf_dr.xelem(42) = 0;
          dHf_dr.xelem(43) = 0;
          dHf_dr.xelem(44) = -4*Zeta2;
          dHf_dr.xelem(45) = 4*Zeta3-4*Zeta1;
          dHf_dr.xelem(46) = 0;
          dHf_dr.xelem(47) = 0;
          dHf_dr.xelem(48) = 0;
          dHf_dr.xelem(49) = 4*Zeta3-4*Zeta1;
          dHf_dr.xelem(50) = 0;
          dHf_dr.xelem(51) = 0;
          dHf_dr.xelem(52) = 0;
          dHf_dr.xelem(53) = 4*Zeta3-4*Zeta1;
     }

     static void VectorInterpMatrixDerS(const ColumnVector& rv, Matrix& dHf_ds) {
          FEM_ASSERT(rv.rows() == 3);
          FEM_ASSERT(dHf_ds.rows() == 3);
          FEM_ASSERT(dHf_ds.columns() == 18);

          const double Zeta1 = rv.xelem(0);
          const double Zeta2 = rv.xelem(1);
          const double Zeta3 = rv.xelem(2);

          dHf_ds.xelem(0) = 0;
          dHf_ds.xelem(1) = 0;
          dHf_ds.xelem(2) = 0;
          dHf_ds.xelem(3) = 0;
          dHf_ds.xelem(4) = 0;
          dHf_ds.xelem(5) = 0;
          dHf_ds.xelem(6) = 0;
          dHf_ds.xelem(7) = 0;
          dHf_ds.xelem(8) = 0;
          dHf_ds.xelem(9) = 4*Zeta2-1;
          dHf_ds.xelem(10) = 0;
          dHf_ds.xelem(11) = 0;
          dHf_ds.xelem(12) = 0;
          dHf_ds.xelem(13) = 4*Zeta2-1;
          dHf_ds.xelem(14) = 0;
          dHf_ds.xelem(15) = 0;
          dHf_ds.xelem(16) = 0;
          dHf_ds.xelem(17) = 4*Zeta2-1;
          dHf_ds.xelem(18) = 1-4*Zeta3;
          dHf_ds.xelem(19) = 0;
          dHf_ds.xelem(20) = 0;
          dHf_ds.xelem(21) = 0;
          dHf_ds.xelem(22) = 1-4*Zeta3;
          dHf_ds.xelem(23) = 0;
          dHf_ds.xelem(24) = 0;
          dHf_ds.xelem(25) = 0;
          dHf_ds.xelem(26) = 1-4*Zeta3;
          dHf_ds.xelem(27) = 4*Zeta1;
          dHf_ds.xelem(28) = 0;
          dHf_ds.xelem(29) = 0;
          dHf_ds.xelem(30) = 0;
          dHf_ds.xelem(31) = 4*Zeta1;
          dHf_ds.xelem(32) = 0;
          dHf_ds.xelem(33) = 0;
          dHf_ds.xelem(34) = 0;
          dHf_ds.xelem(35) = 4*Zeta1;
          dHf_ds.xelem(36) = 4*Zeta3-4*Zeta2;
          dHf_ds.xelem(37) = 0;
          dHf_ds.xelem(38) = 0;
          dHf_ds.xelem(39) = 0;
          dHf_ds.xelem(40) = 4*Zeta3-4*Zeta2;
          dHf_ds.xelem(41) = 0;
          dHf_ds.xelem(42) = 0;
          dHf_ds.xelem(43) = 0;
          dHf_ds.xelem(44) = 4*Zeta3-4*Zeta2;
          dHf_ds.xelem(45) = -4*Zeta1;
          dHf_ds.xelem(46) = 0;
          dHf_ds.xelem(47) = 0;
          dHf_ds.xelem(48) = 0;
          dHf_ds.xelem(49) = -4*Zeta1;
          dHf_ds.xelem(50) = 0;
          dHf_ds.xelem(51) = 0;
          dHf_ds.xelem(52) = 0;
          dHf_ds.xelem(53) = -4*Zeta1;
     }

     static void AllocIntegrationRule(Element::FemMatrixType eMatType) {
          constexpr double tria_area = 0.5; // Factor for triangular area

          switch (eMatType & ~Element::MAT_COLL_PNT_OUTPUT) {
          case Element::VEC_LOAD_LUMPED:
               if (!oIntegLumped.iGetNumEvalPoints()) {
                    oIntegLumped.SetNumEvalPoints(6, 3);

                    constexpr double w1 = 1. / 6.;
                    constexpr double alpha = 1.;
                    constexpr double beta = 0.;

                    for (octave_idx_type i = 0; i < 3; ++i) {
                         oIntegLumped.SetWeight(i, tria_area * w1);

                         for (octave_idx_type j = 0; j < 3; ++j) {
                              oIntegLumped.SetPosition(i, j, i == j ? alpha : beta);
                         }
                    }

                    constexpr double w2 = 1. / 6.;
                    static constexpr double g2[][3] = {{0.5, 0.5, 0.0},
                                                   {0.0, 0.5, 0.5},
                                                   {0.5, 0.0, 0.5}};

                    for (octave_idx_type i = 0; i < 3; ++i) {
                         oIntegLumped.SetWeight(i + 3, tria_area * w2); // Factor 0.5 for triangle

                         for (octave_idx_type j = 0; j < 3; ++j) {
                              oIntegLumped.SetPosition(i + 3, j, g2[i][j]);
                         }
                    }
               }
               break;
          default:
               if (!oIntegConsistent.iGetNumEvalPoints()) {
                    constexpr double g1 = (6. - sqrt(15.)) / 21.;
                    constexpr double g2 = (6. + sqrt(15.)) / 21.;
                    constexpr double g3 = 1. / 3.;
                    constexpr double w1 = (155. - sqrt(15.)) / 1200.;
                    constexpr double w2 = (155. + sqrt(15.)) / 1200.;
                    constexpr double w3 = 9. / 40.;

                    oIntegConsistent.SetNumEvalPoints(7, 3);

                    for (octave_idx_type i = 0; i < 3; ++i) {
                         oIntegConsistent.SetWeight(i, tria_area * w1);

                         for (octave_idx_type j = 0; j < 3; ++j) {
                              oIntegConsistent.SetPosition(i, j, i == j ? 1. - 2. * g1 : g1);
                         }
                    }

                    for (octave_idx_type i = 0; i < 3; ++i) {
                         oIntegConsistent.SetWeight(i + 3, tria_area * w2);

                         for (octave_idx_type j = 0; j < 3; ++j) {
                              oIntegConsistent.SetPosition(i + 3, j, i == j ? 1. - 2. * g2 : g2);
                         }
                    }

                    oIntegConsistent.SetWeight(6, tria_area * w3);

                    for (octave_idx_type j = 0; j < 3; ++j) {
                         oIntegConsistent.SetPosition(6, j, g3);
                    }
               }
          }
     }

     static const IntegrationRule& GetIntegrationRule(Element::FemMatrixType eMatType) {
          const IntegrationRule& oIntegRule = SelectIntegrationRule(eMatType);

          FEM_ASSERT(oIntegRule.iGetNumEvalPoints() > 0);

          return oIntegRule;
     }

private:
     static const IntegrationRule& SelectIntegrationRule(Element::FemMatrixType eMatType) {
          switch (eMatType) {
          case Element::VEC_LOAD_LUMPED:
               return oIntegLumped;
          default:
               return oIntegConsistent;
          }
     }

     static IntegrationRule oIntegLumped, oIntegConsistent;
};

IntegrationRule ShapeTria6::oIntegLumped;
IntegrationRule ShapeTria6::oIntegConsistent;

class ShapeIso4 {
public:
     static constexpr octave_idx_type iGetNumNodes() {
          return 4;
     }

     static constexpr octave_idx_type iGetNumDirections() {
          return 2;
     }

     static constexpr octave_idx_type iGetNumDofNode() {
          return 3;
     }

     static constexpr octave_idx_type iGetNumEqualityConstr() {
          return 0;
     }

     static constexpr octave_idx_type iGetNumInequalityConstr() {
          return 0;
     }

     static constexpr double EqualityConstr(const ColumnVector&) {
          return 0;
     }

     static double InequalityConstr(const ColumnVector& rv) {
          return 0.;
     }

     static void GetElemLimits(ColumnVector& rmin, ColumnVector& rmax) {
          FEM_ASSERT(rmin.rows() == iGetNumDirections());
          FEM_ASSERT(rmax.rows() == rmin.rows());

          for (octave_idx_type i = 0; i < rmin.rows(); ++i) {
               rmin.xelem(i) = -1.;
          }

          for (octave_idx_type i = 0; i < rmax.rows(); ++i) {
               rmax.xelem(i) = 1.;
          }
     }

     static void ScalarInterpMatrix(const ColumnVector& rv, Matrix& HA, octave_idx_type irow) {
          FEM_ASSERT(rv.rows() == 2);
          FEM_ASSERT(HA.rows() > irow);
          FEM_ASSERT(HA.columns() == 4);

          const double r = rv.xelem(0);
          const double s = rv.xelem(1);
          const octave_idx_type nrows = HA.rows();

          HA.xelem(irow + nrows * 0) = ((r+1)*(s+1))/4.0;
          HA.xelem(irow + nrows * 1) = ((1-r)*(s+1))/4.0;
          HA.xelem(irow + nrows * 2) = ((1-r)*(1-s))/4.0;
          HA.xelem(irow + nrows * 3) = ((r+1)*(1-s))/4.0;
     }

     static void VectorInterpMatrix(const ColumnVector& rv, Matrix& Hf) {
          FEM_ASSERT(rv.rows() == 2);
          FEM_ASSERT(Hf.rows() == 3);
          FEM_ASSERT(Hf.columns() == 12);

          const double r = rv.xelem(0);
          const double s = rv.xelem(1);

          Hf.xelem(0) = ((r+1)*(s+1))/4.0E+0;
          Hf.xelem(1) = 0;
          Hf.xelem(2) = 0;
          Hf.xelem(3) = 0;
          Hf.xelem(4) = ((r+1)*(s+1))/4.0E+0;
          Hf.xelem(5) = 0;
          Hf.xelem(6) = 0;
          Hf.xelem(7) = 0;
          Hf.xelem(8) = ((r+1)*(s+1))/4.0E+0;
          Hf.xelem(9) = ((1-r)*(s+1))/4.0E+0;
          Hf.xelem(10) = 0;
          Hf.xelem(11) = 0;
          Hf.xelem(12) = 0;
          Hf.xelem(13) = ((1-r)*(s+1))/4.0E+0;
          Hf.xelem(14) = 0;
          Hf.xelem(15) = 0;
          Hf.xelem(16) = 0;
          Hf.xelem(17) = ((1-r)*(s+1))/4.0E+0;
          Hf.xelem(18) = ((1-r)*(1-s))/4.0E+0;
          Hf.xelem(19) = 0;
          Hf.xelem(20) = 0;
          Hf.xelem(21) = 0;
          Hf.xelem(22) = ((1-r)*(1-s))/4.0E+0;
          Hf.xelem(23) = 0;
          Hf.xelem(24) = 0;
          Hf.xelem(25) = 0;
          Hf.xelem(26) = ((1-r)*(1-s))/4.0E+0;
          Hf.xelem(27) = ((r+1)*(1-s))/4.0E+0;
          Hf.xelem(28) = 0;
          Hf.xelem(29) = 0;
          Hf.xelem(30) = 0;
          Hf.xelem(31) = ((r+1)*(1-s))/4.0E+0;
          Hf.xelem(32) = 0;
          Hf.xelem(33) = 0;
          Hf.xelem(34) = 0;
          Hf.xelem(35) = ((r+1)*(1-s))/4.0E+0;
     }

     static void VectorInterpMatrixDerR(const ColumnVector& rv, Matrix& dHf_dr) {
          FEM_ASSERT(rv.rows() == 2);

          // const double r = rv(0);
          const double s = rv.xelem(1);

          dHf_dr.xelem(0) = (s+1)/4.0E+0;
          dHf_dr.xelem(1) = 0;
          dHf_dr.xelem(2) = 0;
          dHf_dr.xelem(3) = 0;
          dHf_dr.xelem(4) = (s+1)/4.0E+0;
          dHf_dr.xelem(5) = 0;
          dHf_dr.xelem(6) = 0;
          dHf_dr.xelem(7) = 0;
          dHf_dr.xelem(8) = (s+1)/4.0E+0;
          dHf_dr.xelem(9) = -(s+1)/4.0E+0;
          dHf_dr.xelem(10) = 0;
          dHf_dr.xelem(11) = 0;
          dHf_dr.xelem(12) = 0;
          dHf_dr.xelem(13) = -(s+1)/4.0E+0;
          dHf_dr.xelem(14) = 0;
          dHf_dr.xelem(15) = 0;
          dHf_dr.xelem(16) = 0;
          dHf_dr.xelem(17) = -(s+1)/4.0E+0;
          dHf_dr.xelem(18) = -(1-s)/4.0E+0;
          dHf_dr.xelem(19) = 0;
          dHf_dr.xelem(20) = 0;
          dHf_dr.xelem(21) = 0;
          dHf_dr.xelem(22) = -(1-s)/4.0E+0;
          dHf_dr.xelem(23) = 0;
          dHf_dr.xelem(24) = 0;
          dHf_dr.xelem(25) = 0;
          dHf_dr.xelem(26) = -(1-s)/4.0E+0;
          dHf_dr.xelem(27) = (1-s)/4.0E+0;
          dHf_dr.xelem(28) = 0;
          dHf_dr.xelem(29) = 0;
          dHf_dr.xelem(30) = 0;
          dHf_dr.xelem(31) = (1-s)/4.0E+0;
          dHf_dr.xelem(32) = 0;
          dHf_dr.xelem(33) = 0;
          dHf_dr.xelem(34) = 0;
          dHf_dr.xelem(35) = (1-s)/4.0E+0;
     }

     static void VectorInterpMatrixDerS(const ColumnVector& rv, Matrix& dHf_ds) {
          FEM_ASSERT(rv.rows() == 2);

          const double r = rv.xelem(0);

          dHf_ds.xelem(0) = (r+1)/4.0E+0;
          dHf_ds.xelem(1) = 0;
          dHf_ds.xelem(2) = 0;
          dHf_ds.xelem(3) = 0;
          dHf_ds.xelem(4) = (r+1)/4.0E+0;
          dHf_ds.xelem(5) = 0;
          dHf_ds.xelem(6) = 0;
          dHf_ds.xelem(7) = 0;
          dHf_ds.xelem(8) = (r+1)/4.0E+0;
          dHf_ds.xelem(9) = (1-r)/4.0E+0;
          dHf_ds.xelem(10) = 0;
          dHf_ds.xelem(11) = 0;
          dHf_ds.xelem(12) = 0;
          dHf_ds.xelem(13) = (1-r)/4.0E+0;
          dHf_ds.xelem(14) = 0;
          dHf_ds.xelem(15) = 0;
          dHf_ds.xelem(16) = 0;
          dHf_ds.xelem(17) = (1-r)/4.0E+0;
          dHf_ds.xelem(18) = -(1-r)/4.0E+0;
          dHf_ds.xelem(19) = 0;
          dHf_ds.xelem(20) = 0;
          dHf_ds.xelem(21) = 0;
          dHf_ds.xelem(22) = -(1-r)/4.0E+0;
          dHf_ds.xelem(23) = 0;
          dHf_ds.xelem(24) = 0;
          dHf_ds.xelem(25) = 0;
          dHf_ds.xelem(26) = -(1-r)/4.0E+0;
          dHf_ds.xelem(27) = -(r+1)/4.0E+0;
          dHf_ds.xelem(28) = 0;
          dHf_ds.xelem(29) = 0;
          dHf_ds.xelem(30) = 0;
          dHf_ds.xelem(31) = -(r+1)/4.0E+0;
          dHf_ds.xelem(32) = 0;
          dHf_ds.xelem(33) = 0;
          dHf_ds.xelem(34) = 0;
          dHf_ds.xelem(35) = -(r+1)/4.0E+0;
     }

     static void AllocIntegrationRule(Element::FemMatrixType eMatType) {
          static constexpr octave_idx_type N = 2;
          static constexpr double r[2][N] = {{0.577350269189626, -0.577350269189626}, {1., -1.}};
          static constexpr double alpha[2][N] = {{1., 1.}, {1., 1.}};

          const octave_idx_type iIntegRule = SelectIntegrationRule(eMatType);

          if (!rgIntegRule[iIntegRule].iGetNumEvalPoints()) {
               rgIntegRule[iIntegRule].SetNumEvalPoints(N * N, 2);

               octave_idx_type l = 0;

               for (octave_idx_type i = 0; i < N; ++i) {
                    for (octave_idx_type j = 0; j < N; ++j) {
                         rgIntegRule[iIntegRule].SetPosition(l, 0, r[iIntegRule][i]);
                         rgIntegRule[iIntegRule].SetPosition(l, 1, r[iIntegRule][j]);
                         rgIntegRule[iIntegRule].SetWeight(l, alpha[iIntegRule][i] * alpha[iIntegRule][j]);
                         ++l;
                    }
               }
          }
     }

     static const IntegrationRule& GetIntegrationRule(Element::FemMatrixType eMatType) {
          const octave_idx_type iIntegRule = SelectIntegrationRule(eMatType);

          FEM_ASSERT(iIntegRule >= 0);
          FEM_ASSERT(static_cast<size_t>(iIntegRule) < rgIntegRule.size());
          FEM_ASSERT(rgIntegRule[iIntegRule].iGetNumEvalPoints() > 0);

          return rgIntegRule[iIntegRule];
     }

private:
     static octave_idx_type SelectIntegrationRule(Element::FemMatrixType eMatType) {
          switch (eMatType) {
          case Element::VEC_LOAD_LUMPED:
               return 1;
          default:
               return 0;
          }
     }

     static array<IntegrationRule, 2> rgIntegRule;
};

array<IntegrationRule, 2> ShapeIso4::rgIntegRule;

class ShapeQuad8 {
public:
     static constexpr octave_idx_type iGetNumNodes() {
          return 8;
     }

     static constexpr octave_idx_type iGetNumDirections() {
          return 2;
     }

     static constexpr octave_idx_type iGetNumDofNode() {
          return 3;
     }

     static constexpr octave_idx_type iGetNumEqualityConstr() {
          return 0;
     }

     static constexpr octave_idx_type iGetNumInequalityConstr() {
          return 0;
     }

     static constexpr double EqualityConstr(const ColumnVector&) {
          return 0;
     }

     static double InequalityConstr(const ColumnVector& rv) {
          return 0.;
     }

     static void GetElemLimits(ColumnVector& rmin, ColumnVector& rmax) {
          FEM_ASSERT(rmin.rows() == iGetNumDirections());
          FEM_ASSERT(rmax.rows() == rmin.rows());

          for (octave_idx_type i = 0; i < rmin.rows(); ++i) {
               rmin.xelem(i) = -1.;
          }

          for (octave_idx_type i = 0; i < rmax.rows(); ++i) {
               rmax.xelem(i) = 1.;
          }
     }

     static void ScalarInterpMatrix(const ColumnVector& rv, Matrix& HA, octave_idx_type irow) {
          FEM_ASSERT(rv.rows() == 2);
          FEM_ASSERT(HA.rows() > irow);
          FEM_ASSERT(HA.columns() == 8);

          const double r = rv.xelem(0);
          const double s = rv.xelem(1);
          const double r2 = r * r;
          const double s2 = s * s;
          const octave_idx_type nrows = HA.rows();

          HA.xelem(irow + nrows * 0) = ((r+1)*(s+1))/4.0E+0-(((r+1)*(1-s2))/2.0E+0+((1-r2)*(s+1))/2.0E+0)/2.0E+0;
          HA.xelem(irow + nrows * 1) = ((1-r)*(s+1))/4.0E+0-(((1-r)*(1-s2))/2.0E+0+((1-r2)*(s+1))/2.0E+0)/2.0E+0;
          HA.xelem(irow + nrows * 2) = ((1-r)*(1-s))/4.0E+0-(((1-r)*(1-s2))/2.0E+0+((1-r2)*(1-s))/2.0E+0)/2.0E+0;
          HA.xelem(irow + nrows * 3) = ((r+1)*(1-s))/4.0E+0-(((r+1)*(1-s2))/2.0E+0+((1-r2)*(1-s))/2.0E+0)/2.0E+0;
          HA.xelem(irow + nrows * 4) = ((1-r2)*(s+1))/2.0E+0;
          HA.xelem(irow + nrows * 5) = ((1-r)*(1-s2))/2.0E+0;
          HA.xelem(irow + nrows * 6) = ((1-r2)*(1-s))/2.0E+0;
          HA.xelem(irow + nrows * 7) = ((r+1)*(1-s2))/2.0E+0;
     }

     static void VectorInterpMatrix(const ColumnVector& rv, Matrix& Hf) {
          FEM_ASSERT(rv.rows() == 2);

          const double r = rv.xelem(0);
          const double s = rv.xelem(1);
          const double r2 = r * r;
          const double s2 = s * s;

          Hf.xelem(0) = ((r+1)*(s+1))/4.0E+0-(((r+1)*(1-s2))/2.0E+0+((1-r2)*(s+1))/2.0E+0)/2.0E+0;
          Hf.xelem(1) = 0;
          Hf.xelem(2) = 0;
          Hf.xelem(3) = 0;
          Hf.xelem(4) = ((r+1)*(s+1))/4.0E+0-(((r+1)*(1-s2))/2.0E+0+((1-r2)*(s+1))/2.0E+0)/2.0E+0;
          Hf.xelem(5) = 0;
          Hf.xelem(6) = 0;
          Hf.xelem(7) = 0;
          Hf.xelem(8) = ((r+1)*(s+1))/4.0E+0-(((r+1)*(1-s2))/2.0E+0+((1-r2)*(s+1))/2.0E+0)/2.0E+0;
          Hf.xelem(9) = ((1-r)*(s+1))/4.0E+0-(((1-r)*(1-s2))/2.0E+0+((1-r2)*(s+1))/2.0E+0)/2.0E+0;
          Hf.xelem(10) = 0;
          Hf.xelem(11) = 0;
          Hf.xelem(12) = 0;
          Hf.xelem(13) = ((1-r)*(s+1))/4.0E+0-(((1-r)*(1-s2))/2.0E+0+((1-r2)*(s+1))/2.0E+0)/2.0E+0;
          Hf.xelem(14) = 0;
          Hf.xelem(15) = 0;
          Hf.xelem(16) = 0;
          Hf.xelem(17) = ((1-r)*(s+1))/4.0E+0-(((1-r)*(1-s2))/2.0E+0+((1-r2)*(s+1))/2.0E+0)/2.0E+0;
          Hf.xelem(18) = ((1-r)*(1-s))/4.0E+0-(((1-r)*(1-s2))/2.0E+0+((1-r2)*(1-s))/2.0E+0)/2.0E+0;
          Hf.xelem(19) = 0;
          Hf.xelem(20) = 0;
          Hf.xelem(21) = 0;
          Hf.xelem(22) = ((1-r)*(1-s))/4.0E+0-(((1-r)*(1-s2))/2.0E+0+((1-r2)*(1-s))/2.0E+0)/2.0E+0;
          Hf.xelem(23) = 0;
          Hf.xelem(24) = 0;
          Hf.xelem(25) = 0;
          Hf.xelem(26) = ((1-r)*(1-s))/4.0E+0-(((1-r)*(1-s2))/2.0E+0+((1-r2)*(1-s))/2.0E+0)/2.0E+0;
          Hf.xelem(27) = ((r+1)*(1-s))/4.0E+0-(((r+1)*(1-s2))/2.0E+0+((1-r2)*(1-s))/2.0E+0)/2.0E+0;
          Hf.xelem(28) = 0;
          Hf.xelem(29) = 0;
          Hf.xelem(30) = 0;
          Hf.xelem(31) = ((r+1)*(1-s))/4.0E+0-(((r+1)*(1-s2))/2.0E+0+((1-r2)*(1-s))/2.0E+0)/2.0E+0;
          Hf.xelem(32) = 0;
          Hf.xelem(33) = 0;
          Hf.xelem(34) = 0;
          Hf.xelem(35) = ((r+1)*(1-s))/4.0E+0-(((r+1)*(1-s2))/2.0E+0+((1-r2)*(1-s))/2.0E+0)/2.0E+0;
          Hf.xelem(36) = ((1-r2)*(s+1))/2.0E+0;
          Hf.xelem(37) = 0;
          Hf.xelem(38) = 0;
          Hf.xelem(39) = 0;
          Hf.xelem(40) = ((1-r2)*(s+1))/2.0E+0;
          Hf.xelem(41) = 0;
          Hf.xelem(42) = 0;
          Hf.xelem(43) = 0;
          Hf.xelem(44) = ((1-r2)*(s+1))/2.0E+0;
          Hf.xelem(45) = ((1-r)*(1-s2))/2.0E+0;
          Hf.xelem(46) = 0;
          Hf.xelem(47) = 0;
          Hf.xelem(48) = 0;
          Hf.xelem(49) = ((1-r)*(1-s2))/2.0E+0;
          Hf.xelem(50) = 0;
          Hf.xelem(51) = 0;
          Hf.xelem(52) = 0;
          Hf.xelem(53) = ((1-r)*(1-s2))/2.0E+0;
          Hf.xelem(54) = ((1-r2)*(1-s))/2.0E+0;
          Hf.xelem(55) = 0;
          Hf.xelem(56) = 0;
          Hf.xelem(57) = 0;
          Hf.xelem(58) = ((1-r2)*(1-s))/2.0E+0;
          Hf.xelem(59) = 0;
          Hf.xelem(60) = 0;
          Hf.xelem(61) = 0;
          Hf.xelem(62) = ((1-r2)*(1-s))/2.0E+0;
          Hf.xelem(63) = ((r+1)*(1-s2))/2.0E+0;
          Hf.xelem(64) = 0;
          Hf.xelem(65) = 0;
          Hf.xelem(66) = 0;
          Hf.xelem(67) = ((r+1)*(1-s2))/2.0E+0;
          Hf.xelem(68) = 0;
          Hf.xelem(69) = 0;
          Hf.xelem(70) = 0;
          Hf.xelem(71) = ((r+1)*(1-s2))/2.0E+0;
     }

     static void VectorInterpMatrixDerR(const ColumnVector& rv, Matrix& dHf_dr) {
          FEM_ASSERT(rv.rows() == 2);

          const double r = rv.xelem(0);
          const double s = rv.xelem(1);
          const double s2 = s * s;

          dHf_dr.xelem(0) = (s+1)/4.0E+0-((1-s2)/2.0E+0-r*(s+1))/2.0E+0;
          dHf_dr.xelem(1) = 0;
          dHf_dr.xelem(2) = 0;
          dHf_dr.xelem(3) = 0;
          dHf_dr.xelem(4) = (s+1)/4.0E+0-((1-s2)/2.0E+0-r*(s+1))/2.0E+0;
          dHf_dr.xelem(5) = 0;
          dHf_dr.xelem(6) = 0;
          dHf_dr.xelem(7) = 0;
          dHf_dr.xelem(8) = (s+1)/4.0E+0-((1-s2)/2.0E+0-r*(s+1))/2.0E+0;
          dHf_dr.xelem(9) = (-((-(1-s2)/2.0E+0)-r*(s+1))/2.0E+0)-(s+1)/4.0E+0;
          dHf_dr.xelem(10) = 0;
          dHf_dr.xelem(11) = 0;
          dHf_dr.xelem(12) = 0;
          dHf_dr.xelem(13) = (-((-(1-s2)/2.0E+0)-r*(s+1))/2.0E+0)-(s+1)/4.0E+0;
          dHf_dr.xelem(14) = 0;
          dHf_dr.xelem(15) = 0;
          dHf_dr.xelem(16) = 0;
          dHf_dr.xelem(17) = (-((-(1-s2)/2.0E+0)-r*(s+1))/2.0E+0)-(s+1)/4.0E+0;
          dHf_dr.xelem(18) = (-((-(1-s2)/2.0E+0)-r*(1-s))/2.0E+0)-(1-s)/4.0E+0;
          dHf_dr.xelem(19) = 0;
          dHf_dr.xelem(20) = 0;
          dHf_dr.xelem(21) = 0;
          dHf_dr.xelem(22) = (-((-(1-s2)/2.0E+0)-r*(1-s))/2.0E+0)-(1-s)/4.0E+0;
          dHf_dr.xelem(23) = 0;
          dHf_dr.xelem(24) = 0;
          dHf_dr.xelem(25) = 0;
          dHf_dr.xelem(26) = (-((-(1-s2)/2.0E+0)-r*(1-s))/2.0E+0)-(1-s)/4.0E+0;
          dHf_dr.xelem(27) = (1-s)/4.0E+0-((1-s2)/2.0E+0-r*(1-s))/2.0E+0;
          dHf_dr.xelem(28) = 0;
          dHf_dr.xelem(29) = 0;
          dHf_dr.xelem(30) = 0;
          dHf_dr.xelem(31) = (1-s)/4.0E+0-((1-s2)/2.0E+0-r*(1-s))/2.0E+0;
          dHf_dr.xelem(32) = 0;
          dHf_dr.xelem(33) = 0;
          dHf_dr.xelem(34) = 0;
          dHf_dr.xelem(35) = (1-s)/4.0E+0-((1-s2)/2.0E+0-r*(1-s))/2.0E+0;
          dHf_dr.xelem(36) = -r*(s+1);
          dHf_dr.xelem(37) = 0;
          dHf_dr.xelem(38) = 0;
          dHf_dr.xelem(39) = 0;
          dHf_dr.xelem(40) = -r*(s+1);
          dHf_dr.xelem(41) = 0;
          dHf_dr.xelem(42) = 0;
          dHf_dr.xelem(43) = 0;
          dHf_dr.xelem(44) = -r*(s+1);
          dHf_dr.xelem(45) = -(1-s2)/2.0E+0;
          dHf_dr.xelem(46) = 0;
          dHf_dr.xelem(47) = 0;
          dHf_dr.xelem(48) = 0;
          dHf_dr.xelem(49) = -(1-s2)/2.0E+0;
          dHf_dr.xelem(50) = 0;
          dHf_dr.xelem(51) = 0;
          dHf_dr.xelem(52) = 0;
          dHf_dr.xelem(53) = -(1-s2)/2.0E+0;
          dHf_dr.xelem(54) = -r*(1-s);
          dHf_dr.xelem(55) = 0;
          dHf_dr.xelem(56) = 0;
          dHf_dr.xelem(57) = 0;
          dHf_dr.xelem(58) = -r*(1-s);
          dHf_dr.xelem(59) = 0;
          dHf_dr.xelem(60) = 0;
          dHf_dr.xelem(61) = 0;
          dHf_dr.xelem(62) = -r*(1-s);
          dHf_dr.xelem(63) = (1-s2)/2.0E+0;
          dHf_dr.xelem(64) = 0;
          dHf_dr.xelem(65) = 0;
          dHf_dr.xelem(66) = 0;
          dHf_dr.xelem(67) = (1-s2)/2.0E+0;
          dHf_dr.xelem(68) = 0;
          dHf_dr.xelem(69) = 0;
          dHf_dr.xelem(70) = 0;
          dHf_dr.xelem(71) = (1-s2)/2.0E+0;
     }

     static void VectorInterpMatrixDerS(const ColumnVector& rv, Matrix& dHf_ds) {
          FEM_ASSERT(rv.rows() == 2);

          const double r = rv.xelem(0);
          const double s = rv.xelem(1);
          const double r2 = r * r;

          dHf_ds.xelem(0) = (r+1)/4.0E+0-((1-r2)/2.0E+0-(r+1)*s)/2.0E+0;
          dHf_ds.xelem(1) = 0;
          dHf_ds.xelem(2) = 0;
          dHf_ds.xelem(3) = 0;
          dHf_ds.xelem(4) = (r+1)/4.0E+0-((1-r2)/2.0E+0-(r+1)*s)/2.0E+0;
          dHf_ds.xelem(5) = 0;
          dHf_ds.xelem(6) = 0;
          dHf_ds.xelem(7) = 0;
          dHf_ds.xelem(8) = (r+1)/4.0E+0-((1-r2)/2.0E+0-(r+1)*s)/2.0E+0;
          dHf_ds.xelem(9) = (1-r)/4.0E+0-((1-r2)/2.0E+0-(1-r)*s)/2.0E+0;
          dHf_ds.xelem(10) = 0;
          dHf_ds.xelem(11) = 0;
          dHf_ds.xelem(12) = 0;
          dHf_ds.xelem(13) = (1-r)/4.0E+0-((1-r2)/2.0E+0-(1-r)*s)/2.0E+0;
          dHf_ds.xelem(14) = 0;
          dHf_ds.xelem(15) = 0;
          dHf_ds.xelem(16) = 0;
          dHf_ds.xelem(17) = (1-r)/4.0E+0-((1-r2)/2.0E+0-(1-r)*s)/2.0E+0;
          dHf_ds.xelem(18) = (-((-(1-r)*s)-(1-r2)/2.0E+0)/2.0E+0)-(1-r)/4.0E+0;
          dHf_ds.xelem(19) = 0;
          dHf_ds.xelem(20) = 0;
          dHf_ds.xelem(21) = 0;
          dHf_ds.xelem(22) = (-((-(1-r)*s)-(1-r2)/2.0E+0)/2.0E+0)-(1-r)/4.0E+0;
          dHf_ds.xelem(23) = 0;
          dHf_ds.xelem(24) = 0;
          dHf_ds.xelem(25) = 0;
          dHf_ds.xelem(26) = (-((-(1-r)*s)-(1-r2)/2.0E+0)/2.0E+0)-(1-r)/4.0E+0;
          dHf_ds.xelem(27) = (-((-(r+1)*s)-(1-r2)/2.0E+0)/2.0E+0)-(r+1)/4.0E+0;
          dHf_ds.xelem(28) = 0;
          dHf_ds.xelem(29) = 0;
          dHf_ds.xelem(30) = 0;
          dHf_ds.xelem(31) = (-((-(r+1)*s)-(1-r2)/2.0E+0)/2.0E+0)-(r+1)/4.0E+0;
          dHf_ds.xelem(32) = 0;
          dHf_ds.xelem(33) = 0;
          dHf_ds.xelem(34) = 0;
          dHf_ds.xelem(35) = (-((-(r+1)*s)-(1-r2)/2.0E+0)/2.0E+0)-(r+1)/4.0E+0;
          dHf_ds.xelem(36) = (1-r2)/2.0E+0;
          dHf_ds.xelem(37) = 0;
          dHf_ds.xelem(38) = 0;
          dHf_ds.xelem(39) = 0;
          dHf_ds.xelem(40) = (1-r2)/2.0E+0;
          dHf_ds.xelem(41) = 0;
          dHf_ds.xelem(42) = 0;
          dHf_ds.xelem(43) = 0;
          dHf_ds.xelem(44) = (1-r2)/2.0E+0;
          dHf_ds.xelem(45) = -(1-r)*s;
          dHf_ds.xelem(46) = 0;
          dHf_ds.xelem(47) = 0;
          dHf_ds.xelem(48) = 0;
          dHf_ds.xelem(49) = -(1-r)*s;
          dHf_ds.xelem(50) = 0;
          dHf_ds.xelem(51) = 0;
          dHf_ds.xelem(52) = 0;
          dHf_ds.xelem(53) = -(1-r)*s;
          dHf_ds.xelem(54) = -(1-r2)/2.0E+0;
          dHf_ds.xelem(55) = 0;
          dHf_ds.xelem(56) = 0;
          dHf_ds.xelem(57) = 0;
          dHf_ds.xelem(58) = -(1-r2)/2.0E+0;
          dHf_ds.xelem(59) = 0;
          dHf_ds.xelem(60) = 0;
          dHf_ds.xelem(61) = 0;
          dHf_ds.xelem(62) = -(1-r2)/2.0E+0;
          dHf_ds.xelem(63) = -(r+1)*s;
          dHf_ds.xelem(64) = 0;
          dHf_ds.xelem(65) = 0;
          dHf_ds.xelem(66) = 0;
          dHf_ds.xelem(67) = -(r+1)*s;
          dHf_ds.xelem(68) = 0;
          dHf_ds.xelem(69) = 0;
          dHf_ds.xelem(70) = 0;
          dHf_ds.xelem(71) = -(r+1)*s;
     }

     static void AllocIntegrationRule(Element::FemMatrixType eMatType) {
          static constexpr octave_idx_type N = 3;
          static constexpr double r[2][N] = {{0.774596669241483, 0., -0.774596669241483}, {1., 0., -1.}};
          static constexpr double alpha[2][N] = {{0.555555555555556, 0.888888888888889, 0.555555555555556}, {2./3., 2./3., 2./3.}};

          const octave_idx_type iIntegRule = SelectIntegrationRule(eMatType);

          if (!rgIntegRule[iIntegRule].iGetNumEvalPoints()) {
               rgIntegRule[iIntegRule].SetNumEvalPoints(N * N, 2);

               octave_idx_type l = 0;

               for (octave_idx_type i = 0; i < N; ++i) {
                    for (octave_idx_type j = 0; j < N; ++j) {
                         rgIntegRule[iIntegRule].SetPosition(l, 0, r[iIntegRule][i]);
                         rgIntegRule[iIntegRule].SetPosition(l, 1, r[iIntegRule][j]);
                         rgIntegRule[iIntegRule].SetWeight(l, alpha[iIntegRule][i] * alpha[iIntegRule][j]);
                         ++l;
                    }
               }
          }
     }

     static const IntegrationRule& GetIntegrationRule(Element::FemMatrixType eMatType) {
          const octave_idx_type iIntegRule = SelectIntegrationRule(eMatType);

          FEM_ASSERT(iIntegRule >= 0);
          FEM_ASSERT(static_cast<size_t>(iIntegRule) < rgIntegRule.size());
          FEM_ASSERT(rgIntegRule[iIntegRule].iGetNumEvalPoints() > 0);

          return rgIntegRule[iIntegRule];
     }

private:
    static octave_idx_type SelectIntegrationRule(Element::FemMatrixType eMatType) {
          switch (eMatType) {
          case Element::VEC_LOAD_LUMPED:
               return 1;
          default:
               return 0;
          }
     }

     static array<IntegrationRule, 2> rgIntegRule;
};

array<IntegrationRule, 2> ShapeQuad8::rgIntegRule;

class ShapeQuad8r {
public:
     static constexpr octave_idx_type iGetNumNodes() {
          return 8;
     }

     static constexpr octave_idx_type iGetNumDirections() {
          return 2;
     }

     static constexpr octave_idx_type iGetNumDofNode() {
          return 3;
     }

     static constexpr octave_idx_type iGetNumEqualityConstr() {
          return 0;
     }

     static constexpr octave_idx_type iGetNumInequalityConstr() {
          return 0;
     }

     static constexpr double EqualityConstr(const ColumnVector&) {
          return 0;
     }

     static double InequalityConstr(const ColumnVector& rv) {
          return 0.;
     }

     static void GetElemLimits(ColumnVector& rmin, ColumnVector& rmax) {
          FEM_ASSERT(rmin.rows() == iGetNumDirections());
          FEM_ASSERT(rmax.rows() == rmin.rows());

          for (octave_idx_type i = 0; i < rmin.rows(); ++i) {
               rmin.xelem(i) = -1.;
          }

          for (octave_idx_type i = 0; i < rmax.rows(); ++i) {
               rmax.xelem(i) = 1.;
          }
     }

     static void ScalarInterpMatrix(const ColumnVector& rv, Matrix& HA, octave_idx_type irow) {
          FEM_ASSERT(rv.rows() == 2);
          FEM_ASSERT(HA.rows() > irow);
          FEM_ASSERT(HA.columns() == 8);

          const double r1 = rv.xelem(0);
          const double r2 = rv.xelem(1);
          const octave_idx_type nrows = HA.rows();

          HA.xelem(irow + nrows * 0) = 2.5E-1*(r1-1.0E+0)*(1.0E+0-r2)*(r2+r1+1.0E+0);
          HA.xelem(irow + nrows * 1) = 2.5E-1*((-r1)-1.0E+0)*(1.0E+0-r2)*(r2-r1+1.0E+0);
          HA.xelem(irow + nrows * 2) = 2.5E-1*((-r1)-1.0E+0)*((-r2)-r1+1.0E+0)*(r2+1.0E+0);
          HA.xelem(irow + nrows * 3) = 2.5E-1*(r1-1.0E+0)*((-r2)+r1+1.0E+0)*(r2+1.0E+0);
          HA.xelem(irow + nrows * 4) = 5.0E-1*(1.0E+0-r1)*(r1+1.0E+0)*(1.0E+0-r2);
          HA.xelem(irow + nrows * 5) = 5.0E-1*(r1+1.0E+0)*(1.0E+0-r2)*(r2+1.0E+0);
          HA.xelem(irow + nrows * 6) = 5.0E-1*(1.0E+0-r1)*(r1+1.0E+0)*(r2+1.0E+0);
          HA.xelem(irow + nrows * 7) = 5.0E-1*(1.0E+0-r1)*(1.0E+0-r2)*(r2+1.0E+0);
     }

     static void VectorInterpMatrix(const ColumnVector& rv, Matrix& Hf) {
          FEM_ASSERT(rv.rows() == 2);

          const double r1 = rv.xelem(0);
          const double r2 = rv.xelem(1);

          Hf.xelem(0) = 2.5E-1*(r1-1.0E+0)*(1.0E+0-r2)*(r2+r1+1.0E+0);
          Hf.xelem(1) = 0;
          Hf.xelem(2) = 0;
          Hf.xelem(3) = 0;
          Hf.xelem(4) = 2.5E-1*(r1-1.0E+0)*(1.0E+0-r2)*(r2+r1+1.0E+0);
          Hf.xelem(5) = 0;
          Hf.xelem(6) = 0;
          Hf.xelem(7) = 0;
          Hf.xelem(8) = 2.5E-1*(r1-1.0E+0)*(1.0E+0-r2)*(r2+r1+1.0E+0);
          Hf.xelem(9) = 2.5E-1*((-r1)-1.0E+0)*(1.0E+0-r2)*(r2-r1+1.0E+0);
          Hf.xelem(10) = 0;
          Hf.xelem(11) = 0;
          Hf.xelem(12) = 0;
          Hf.xelem(13) = 2.5E-1*((-r1)-1.0E+0)*(1.0E+0-r2)*(r2-r1+1.0E+0);
          Hf.xelem(14) = 0;
          Hf.xelem(15) = 0;
          Hf.xelem(16) = 0;
          Hf.xelem(17) = 2.5E-1*((-r1)-1.0E+0)*(1.0E+0-r2)*(r2-r1+1.0E+0);
          Hf.xelem(18) = 2.5E-1*((-r1)-1.0E+0)*((-r2)-r1+1.0E+0)*(r2+1.0E+0);
          Hf.xelem(19) = 0;
          Hf.xelem(20) = 0;
          Hf.xelem(21) = 0;
          Hf.xelem(22) = 2.5E-1*((-r1)-1.0E+0)*((-r2)-r1+1.0E+0)*(r2+1.0E+0);
          Hf.xelem(23) = 0;
          Hf.xelem(24) = 0;
          Hf.xelem(25) = 0;
          Hf.xelem(26) = 2.5E-1*((-r1)-1.0E+0)*((-r2)-r1+1.0E+0)*(r2+1.0E+0);
          Hf.xelem(27) = 2.5E-1*(r1-1.0E+0)*((-r2)+r1+1.0E+0)*(r2+1.0E+0);
          Hf.xelem(28) = 0;
          Hf.xelem(29) = 0;
          Hf.xelem(30) = 0;
          Hf.xelem(31) = 2.5E-1*(r1-1.0E+0)*((-r2)+r1+1.0E+0)*(r2+1.0E+0);
          Hf.xelem(32) = 0;
          Hf.xelem(33) = 0;
          Hf.xelem(34) = 0;
          Hf.xelem(35) = 2.5E-1*(r1-1.0E+0)*((-r2)+r1+1.0E+0)*(r2+1.0E+0);
          Hf.xelem(36) = 5.0E-1*(1.0E+0-r1)*(r1+1.0E+0)*(1.0E+0-r2);
          Hf.xelem(37) = 0;
          Hf.xelem(38) = 0;
          Hf.xelem(39) = 0;
          Hf.xelem(40) = 5.0E-1*(1.0E+0-r1)*(r1+1.0E+0)*(1.0E+0-r2);
          Hf.xelem(41) = 0;
          Hf.xelem(42) = 0;
          Hf.xelem(43) = 0;
          Hf.xelem(44) = 5.0E-1*(1.0E+0-r1)*(r1+1.0E+0)*(1.0E+0-r2);
          Hf.xelem(45) = 5.0E-1*(r1+1.0E+0)*(1.0E+0-r2)*(r2+1.0E+0);
          Hf.xelem(46) = 0;
          Hf.xelem(47) = 0;
          Hf.xelem(48) = 0;
          Hf.xelem(49) = 5.0E-1*(r1+1.0E+0)*(1.0E+0-r2)*(r2+1.0E+0);
          Hf.xelem(50) = 0;
          Hf.xelem(51) = 0;
          Hf.xelem(52) = 0;
          Hf.xelem(53) = 5.0E-1*(r1+1.0E+0)*(1.0E+0-r2)*(r2+1.0E+0);
          Hf.xelem(54) = 5.0E-1*(1.0E+0-r1)*(r1+1.0E+0)*(r2+1.0E+0);
          Hf.xelem(55) = 0;
          Hf.xelem(56) = 0;
          Hf.xelem(57) = 0;
          Hf.xelem(58) = 5.0E-1*(1.0E+0-r1)*(r1+1.0E+0)*(r2+1.0E+0);
          Hf.xelem(59) = 0;
          Hf.xelem(60) = 0;
          Hf.xelem(61) = 0;
          Hf.xelem(62) = 5.0E-1*(1.0E+0-r1)*(r1+1.0E+0)*(r2+1.0E+0);
          Hf.xelem(63) = 5.0E-1*(1.0E+0-r1)*(1.0E+0-r2)*(r2+1.0E+0);
          Hf.xelem(64) = 0;
          Hf.xelem(65) = 0;
          Hf.xelem(66) = 0;
          Hf.xelem(67) = 5.0E-1*(1.0E+0-r1)*(1.0E+0-r2)*(r2+1.0E+0);
          Hf.xelem(68) = 0;
          Hf.xelem(69) = 0;
          Hf.xelem(70) = 0;
          Hf.xelem(71) = 5.0E-1*(1.0E+0-r1)*(1.0E+0-r2)*(r2+1.0E+0);



     }

     static void VectorInterpMatrixDerR(const ColumnVector& rv, Matrix& dHf_dr) {
          FEM_ASSERT(rv.rows() == 2);

          const double r1 = rv.xelem(0);
          const double r2 = rv.xelem(1);

          dHf_dr.xelem(0) = 2.5E-1*(1.0E+0-r2)*(r2+r1+1.0E+0)+2.5E-1*(r1-1.0E+0)*(1.0E+0-r2);
          dHf_dr.xelem(1) = 0;
          dHf_dr.xelem(2) = 0;
          dHf_dr.xelem(3) = 0;
          dHf_dr.xelem(4) = 2.5E-1*(1.0E+0-r2)*(r2+r1+1.0E+0)+2.5E-1*(r1-1.0E+0)*(1.0E+0-r2);
          dHf_dr.xelem(5) = 0;
          dHf_dr.xelem(6) = 0;
          dHf_dr.xelem(7) = 0;
          dHf_dr.xelem(8) = 2.5E-1*(1.0E+0-r2)*(r2+r1+1.0E+0)+2.5E-1*(r1-1.0E+0)*(1.0E+0-r2);
          dHf_dr.xelem(9) = (-2.5E-1*(1.0E+0-r2)*(r2-r1+1.0E+0))-2.5E-1*((-r1)-1.0E+0)*(1.0E+0-r2);
          dHf_dr.xelem(10) = 0;
          dHf_dr.xelem(11) = 0;
          dHf_dr.xelem(12) = 0;
          dHf_dr.xelem(13) = (-2.5E-1*(1.0E+0-r2)*(r2-r1+1.0E+0))-2.5E-1*((-r1)-1.0E+0)*(1.0E+0-r2);
          dHf_dr.xelem(14) = 0;
          dHf_dr.xelem(15) = 0;
          dHf_dr.xelem(16) = 0;
          dHf_dr.xelem(17) = (-2.5E-1*(1.0E+0-r2)*(r2-r1+1.0E+0))-2.5E-1*((-r1)-1.0E+0)*(1.0E+0-r2);
          dHf_dr.xelem(18) = (-2.5E-1*((-r2)-r1+1.0E+0)*(r2+1.0E+0))-2.5E-1*((-r1)-1.0E+0)*(r2+1.0E+0);
          dHf_dr.xelem(19) = 0;
          dHf_dr.xelem(20) = 0;
          dHf_dr.xelem(21) = 0;
          dHf_dr.xelem(22) = (-2.5E-1*((-r2)-r1+1.0E+0)*(r2+1.0E+0))-2.5E-1*((-r1)-1.0E+0)*(r2+1.0E+0);
          dHf_dr.xelem(23) = 0;
          dHf_dr.xelem(24) = 0;
          dHf_dr.xelem(25) = 0;
          dHf_dr.xelem(26) = (-2.5E-1*((-r2)-r1+1.0E+0)*(r2+1.0E+0))-2.5E-1*((-r1)-1.0E+0)*(r2+1.0E+0);
          dHf_dr.xelem(27) = 2.5E-1*((-r2)+r1+1.0E+0)*(r2+1.0E+0)+2.5E-1*(r1-1.0E+0)*(r2+1.0E+0);
          dHf_dr.xelem(28) = 0;
          dHf_dr.xelem(29) = 0;
          dHf_dr.xelem(30) = 0;
          dHf_dr.xelem(31) = 2.5E-1*((-r2)+r1+1.0E+0)*(r2+1.0E+0)+2.5E-1*(r1-1.0E+0)*(r2+1.0E+0);
          dHf_dr.xelem(32) = 0;
          dHf_dr.xelem(33) = 0;
          dHf_dr.xelem(34) = 0;
          dHf_dr.xelem(35) = 2.5E-1*((-r2)+r1+1.0E+0)*(r2+1.0E+0)+2.5E-1*(r1-1.0E+0)*(r2+1.0E+0);
          dHf_dr.xelem(36) = 5.0E-1*(1.0E+0-r1)*(1.0E+0-r2)-5.0E-1*(r1+1.0E+0)*(1.0E+0-r2);
          dHf_dr.xelem(37) = 0;
          dHf_dr.xelem(38) = 0;
          dHf_dr.xelem(39) = 0;
          dHf_dr.xelem(40) = 5.0E-1*(1.0E+0-r1)*(1.0E+0-r2)-5.0E-1*(r1+1.0E+0)*(1.0E+0-r2);
          dHf_dr.xelem(41) = 0;
          dHf_dr.xelem(42) = 0;
          dHf_dr.xelem(43) = 0;
          dHf_dr.xelem(44) = 5.0E-1*(1.0E+0-r1)*(1.0E+0-r2)-5.0E-1*(r1+1.0E+0)*(1.0E+0-r2);
          dHf_dr.xelem(45) = 5.0E-1*(1.0E+0-r2)*(r2+1.0E+0);
          dHf_dr.xelem(46) = 0;
          dHf_dr.xelem(47) = 0;
          dHf_dr.xelem(48) = 0;
          dHf_dr.xelem(49) = 5.0E-1*(1.0E+0-r2)*(r2+1.0E+0);
          dHf_dr.xelem(50) = 0;
          dHf_dr.xelem(51) = 0;
          dHf_dr.xelem(52) = 0;
          dHf_dr.xelem(53) = 5.0E-1*(1.0E+0-r2)*(r2+1.0E+0);
          dHf_dr.xelem(54) = 5.0E-1*(1.0E+0-r1)*(r2+1.0E+0)-5.0E-1*(r1+1.0E+0)*(r2+1.0E+0);
          dHf_dr.xelem(55) = 0;
          dHf_dr.xelem(56) = 0;
          dHf_dr.xelem(57) = 0;
          dHf_dr.xelem(58) = 5.0E-1*(1.0E+0-r1)*(r2+1.0E+0)-5.0E-1*(r1+1.0E+0)*(r2+1.0E+0);
          dHf_dr.xelem(59) = 0;
          dHf_dr.xelem(60) = 0;
          dHf_dr.xelem(61) = 0;
          dHf_dr.xelem(62) = 5.0E-1*(1.0E+0-r1)*(r2+1.0E+0)-5.0E-1*(r1+1.0E+0)*(r2+1.0E+0);
          dHf_dr.xelem(63) = -5.0E-1*(1.0E+0-r2)*(r2+1.0E+0);
          dHf_dr.xelem(64) = 0;
          dHf_dr.xelem(65) = 0;
          dHf_dr.xelem(66) = 0;
          dHf_dr.xelem(67) = -5.0E-1*(1.0E+0-r2)*(r2+1.0E+0);
          dHf_dr.xelem(68) = 0;
          dHf_dr.xelem(69) = 0;
          dHf_dr.xelem(70) = 0;
          dHf_dr.xelem(71) = -5.0E-1*(1.0E+0-r2)*(r2+1.0E+0);
     }

     static void VectorInterpMatrixDerS(const ColumnVector& rv, Matrix& dHf_ds) {
          FEM_ASSERT(rv.rows() == 2);

          const double r1 = rv.xelem(0);
          const double r2 = rv.xelem(1);

          dHf_ds.xelem(0) = 2.5E-1*(r1-1.0E+0)*(1.0E+0-r2)-2.5E-1*(r1-1.0E+0)*(r2+r1+1.0E+0);
          dHf_ds.xelem(1) = 0;
          dHf_ds.xelem(2) = 0;
          dHf_ds.xelem(3) = 0;
          dHf_ds.xelem(4) = 2.5E-1*(r1-1.0E+0)*(1.0E+0-r2)-2.5E-1*(r1-1.0E+0)*(r2+r1+1.0E+0);
          dHf_ds.xelem(5) = 0;
          dHf_ds.xelem(6) = 0;
          dHf_ds.xelem(7) = 0;
          dHf_ds.xelem(8) = 2.5E-1*(r1-1.0E+0)*(1.0E+0-r2)-2.5E-1*(r1-1.0E+0)*(r2+r1+1.0E+0);
          dHf_ds.xelem(9) = 2.5E-1*((-r1)-1.0E+0)*(1.0E+0-r2)-2.5E-1*((-r1)-1.0E+0)*(r2-r1+1.0E+0);
          dHf_ds.xelem(10) = 0;
          dHf_ds.xelem(11) = 0;
          dHf_ds.xelem(12) = 0;
          dHf_ds.xelem(13) = 2.5E-1*((-r1)-1.0E+0)*(1.0E+0-r2)-2.5E-1*((-r1)-1.0E+0)*(r2-r1+1.0E+0);
          dHf_ds.xelem(14) = 0;
          dHf_ds.xelem(15) = 0;
          dHf_ds.xelem(16) = 0;
          dHf_ds.xelem(17) = 2.5E-1*((-r1)-1.0E+0)*(1.0E+0-r2)-2.5E-1*((-r1)-1.0E+0)*(r2-r1+1.0E+0);
          dHf_ds.xelem(18) = 2.5E-1*((-r1)-1.0E+0)*((-r2)-r1+1.0E+0)-2.5E-1*((-r1)-1.0E+0)*(r2+1.0E+0);
          dHf_ds.xelem(19) = 0;
          dHf_ds.xelem(20) = 0;
          dHf_ds.xelem(21) = 0;
          dHf_ds.xelem(22) = 2.5E-1*((-r1)-1.0E+0)*((-r2)-r1+1.0E+0)-2.5E-1*((-r1)-1.0E+0)*(r2+1.0E+0);
          dHf_ds.xelem(23) = 0;
          dHf_ds.xelem(24) = 0;
          dHf_ds.xelem(25) = 0;
          dHf_ds.xelem(26) = 2.5E-1*((-r1)-1.0E+0)*((-r2)-r1+1.0E+0)-2.5E-1*((-r1)-1.0E+0)*(r2+1.0E+0);
          dHf_ds.xelem(27) = 2.5E-1*(r1-1.0E+0)*((-r2)+r1+1.0E+0)-2.5E-1*(r1-1.0E+0)*(r2+1.0E+0);
          dHf_ds.xelem(28) = 0;
          dHf_ds.xelem(29) = 0;
          dHf_ds.xelem(30) = 0;
          dHf_ds.xelem(31) = 2.5E-1*(r1-1.0E+0)*((-r2)+r1+1.0E+0)-2.5E-1*(r1-1.0E+0)*(r2+1.0E+0);
          dHf_ds.xelem(32) = 0;
          dHf_ds.xelem(33) = 0;
          dHf_ds.xelem(34) = 0;
          dHf_ds.xelem(35) = 2.5E-1*(r1-1.0E+0)*((-r2)+r1+1.0E+0)-2.5E-1*(r1-1.0E+0)*(r2+1.0E+0);
          dHf_ds.xelem(36) = -5.0E-1*(1.0E+0-r1)*(r1+1.0E+0);
          dHf_ds.xelem(37) = 0;
          dHf_ds.xelem(38) = 0;
          dHf_ds.xelem(39) = 0;
          dHf_ds.xelem(40) = -5.0E-1*(1.0E+0-r1)*(r1+1.0E+0);
          dHf_ds.xelem(41) = 0;
          dHf_ds.xelem(42) = 0;
          dHf_ds.xelem(43) = 0;
          dHf_ds.xelem(44) = -5.0E-1*(1.0E+0-r1)*(r1+1.0E+0);
          dHf_ds.xelem(45) = 5.0E-1*(r1+1.0E+0)*(1.0E+0-r2)-5.0E-1*(r1+1.0E+0)*(r2+1.0E+0);
          dHf_ds.xelem(46) = 0;
          dHf_ds.xelem(47) = 0;
          dHf_ds.xelem(48) = 0;
          dHf_ds.xelem(49) = 5.0E-1*(r1+1.0E+0)*(1.0E+0-r2)-5.0E-1*(r1+1.0E+0)*(r2+1.0E+0);
          dHf_ds.xelem(50) = 0;
          dHf_ds.xelem(51) = 0;
          dHf_ds.xelem(52) = 0;
          dHf_ds.xelem(53) = 5.0E-1*(r1+1.0E+0)*(1.0E+0-r2)-5.0E-1*(r1+1.0E+0)*(r2+1.0E+0);
          dHf_ds.xelem(54) = 5.0E-1*(1.0E+0-r1)*(r1+1.0E+0);
          dHf_ds.xelem(55) = 0;
          dHf_ds.xelem(56) = 0;
          dHf_ds.xelem(57) = 0;
          dHf_ds.xelem(58) = 5.0E-1*(1.0E+0-r1)*(r1+1.0E+0);
          dHf_ds.xelem(59) = 0;
          dHf_ds.xelem(60) = 0;
          dHf_ds.xelem(61) = 0;
          dHf_ds.xelem(62) = 5.0E-1*(1.0E+0-r1)*(r1+1.0E+0);
          dHf_ds.xelem(63) = 5.0E-1*(1.0E+0-r1)*(1.0E+0-r2)-5.0E-1*(1.0E+0-r1)*(r2+1.0E+0);
          dHf_ds.xelem(64) = 0;
          dHf_ds.xelem(65) = 0;
          dHf_ds.xelem(66) = 0;
          dHf_ds.xelem(67) = 5.0E-1*(1.0E+0-r1)*(1.0E+0-r2)-5.0E-1*(1.0E+0-r1)*(r2+1.0E+0);
          dHf_ds.xelem(68) = 0;
          dHf_ds.xelem(69) = 0;
          dHf_ds.xelem(70) = 0;
          dHf_ds.xelem(71) = 5.0E-1*(1.0E+0-r1)*(1.0E+0-r2)-5.0E-1*(1.0E+0-r1)*(r2+1.0E+0);
     }

     static void AllocIntegrationRule(Element::FemMatrixType eMatType) {
          static constexpr octave_idx_type N = 3;
          static constexpr double r[2][N] = {{0.774596669241483, 0., -0.774596669241483}, {1., 0., -1.}};
          static constexpr double alpha[2][N] = {{0.555555555555556, 0.888888888888889, 0.555555555555556}, {2./3., 2./3., 2./3.}};

          const octave_idx_type iIntegRule = SelectIntegrationRule(eMatType);

          if (!rgIntegRule[iIntegRule].iGetNumEvalPoints()) {
               rgIntegRule[iIntegRule].SetNumEvalPoints(N * N, 2);

               octave_idx_type l = 0;

               for (octave_idx_type i = 0; i < N; ++i) {
                    for (octave_idx_type j = 0; j < N; ++j) {
                         rgIntegRule[iIntegRule].SetPosition(l, 0, r[iIntegRule][i]);
                         rgIntegRule[iIntegRule].SetPosition(l, 1, r[iIntegRule][j]);
                         rgIntegRule[iIntegRule].SetWeight(l, alpha[iIntegRule][i] * alpha[iIntegRule][j]);
                         ++l;
                    }
               }
          }
     }

     static const IntegrationRule& GetIntegrationRule(Element::FemMatrixType eMatType) {
          const octave_idx_type iIntegRule = SelectIntegrationRule(eMatType);

          FEM_ASSERT(iIntegRule >= 0);
          FEM_ASSERT(static_cast<size_t>(iIntegRule) < rgIntegRule.size());
          FEM_ASSERT(rgIntegRule[iIntegRule].iGetNumEvalPoints() > 0);

          return rgIntegRule[iIntegRule];
     }

private:
     static octave_idx_type SelectIntegrationRule(Element::FemMatrixType eMatType) {
          switch (eMatType) {
          case Element::VEC_LOAD_LUMPED:
               return 1;
          default:
               return 0;
          }
     }

     static array<IntegrationRule, 2> rgIntegRule;
};

array<IntegrationRule, 2> ShapeQuad8r::rgIntegRule;

class ShapeQuad9 {
public:
     static constexpr octave_idx_type iGetNumNodes() {
          return 9;
     }

     static constexpr octave_idx_type iGetNumDirections() {
          return 2;
     }

     static constexpr octave_idx_type iGetNumDofNode() {
          return 3;
     }

     static constexpr octave_idx_type iGetNumEqualityConstr() {
          return 0;
     }

     static constexpr octave_idx_type iGetNumInequalityConstr() {
          return 0;
     }

     static constexpr double EqualityConstr(const ColumnVector&) {
          return 0;
     }

     static double InequalityConstr(const ColumnVector& rv) {
          return 0.;
     }

     static void GetElemLimits(ColumnVector& rmin, ColumnVector& rmax) {
          FEM_ASSERT(rmin.rows() == iGetNumDirections());
          FEM_ASSERT(rmax.rows() == rmin.rows());

          for (octave_idx_type i = 0; i < rmin.rows(); ++i) {
               rmin.xelem(i) = -1.;
          }

          for (octave_idx_type i = 0; i < rmax.rows(); ++i) {
               rmax.xelem(i) = 1.;
          }
     }

     static void ScalarInterpMatrix(const ColumnVector& rv, Matrix& HA, octave_idx_type irow) {
          FEM_ASSERT(rv.rows() == 2);
          FEM_ASSERT(HA.rows() > irow);
          FEM_ASSERT(HA.columns() == 9);

          const double r = rv.xelem(0);
          const double s = rv.xelem(1);
          const double r2 = r * r;
          const double s2 = s * s;
          const octave_idx_type nrows = HA.rows();

          HA.xelem(irow + nrows * 0) = ((r-1)*r*(s-1)*s)/4.0E+0;
          HA.xelem(irow + nrows * 1) = (r*(r+1)*(s-1)*s)/4.0E+0;
          HA.xelem(irow + nrows * 2) = (r*(r+1)*s*(s+1))/4.0E+0;
          HA.xelem(irow + nrows * 3) = ((r-1)*r*s*(s+1))/4.0E+0;
          HA.xelem(irow + nrows * 4) = ((1-r2)*(s-1)*s)/2.0E+0;
          HA.xelem(irow + nrows * 5) = (r*(r+1)*(1-s2))/2.0E+0;
          HA.xelem(irow + nrows * 6) = ((1-r2)*s*(s+1))/2.0E+0;
          HA.xelem(irow + nrows * 7) = ((r-1)*r*(1-s2))/2.0E+0;
          HA.xelem(irow + nrows * 8) = (1-r2)*(1-s2);
     }

     static void VectorInterpMatrix(const ColumnVector& rv, Matrix& Hf) {
          FEM_ASSERT(rv.rows() == 2);

          const double r = rv.xelem(0);
          const double s = rv.xelem(1);
          const double r2 = r * r;
          const double s2 = s * s;


          Hf.xelem(0) = ((r-1)*r*(s-1)*s)/4.0E+0;
          Hf.xelem(1) = 0;
          Hf.xelem(2) = 0;
          Hf.xelem(3) = 0;
          Hf.xelem(4) = ((r-1)*r*(s-1)*s)/4.0E+0;
          Hf.xelem(5) = 0;
          Hf.xelem(6) = 0;
          Hf.xelem(7) = 0;
          Hf.xelem(8) = ((r-1)*r*(s-1)*s)/4.0E+0;
          Hf.xelem(9) = (r*(r+1)*(s-1)*s)/4.0E+0;
          Hf.xelem(10) = 0;
          Hf.xelem(11) = 0;
          Hf.xelem(12) = 0;
          Hf.xelem(13) = (r*(r+1)*(s-1)*s)/4.0E+0;
          Hf.xelem(14) = 0;
          Hf.xelem(15) = 0;
          Hf.xelem(16) = 0;
          Hf.xelem(17) = (r*(r+1)*(s-1)*s)/4.0E+0;
          Hf.xelem(18) = (r*(r+1)*s*(s+1))/4.0E+0;
          Hf.xelem(19) = 0;
          Hf.xelem(20) = 0;
          Hf.xelem(21) = 0;
          Hf.xelem(22) = (r*(r+1)*s*(s+1))/4.0E+0;
          Hf.xelem(23) = 0;
          Hf.xelem(24) = 0;
          Hf.xelem(25) = 0;
          Hf.xelem(26) = (r*(r+1)*s*(s+1))/4.0E+0;
          Hf.xelem(27) = ((r-1)*r*s*(s+1))/4.0E+0;
          Hf.xelem(28) = 0;
          Hf.xelem(29) = 0;
          Hf.xelem(30) = 0;
          Hf.xelem(31) = ((r-1)*r*s*(s+1))/4.0E+0;
          Hf.xelem(32) = 0;
          Hf.xelem(33) = 0;
          Hf.xelem(34) = 0;
          Hf.xelem(35) = ((r-1)*r*s*(s+1))/4.0E+0;
          Hf.xelem(36) = ((1-r2)*(s-1)*s)/2.0E+0;
          Hf.xelem(37) = 0;
          Hf.xelem(38) = 0;
          Hf.xelem(39) = 0;
          Hf.xelem(40) = ((1-r2)*(s-1)*s)/2.0E+0;
          Hf.xelem(41) = 0;
          Hf.xelem(42) = 0;
          Hf.xelem(43) = 0;
          Hf.xelem(44) = ((1-r2)*(s-1)*s)/2.0E+0;
          Hf.xelem(45) = (r*(r+1)*(1-s2))/2.0E+0;
          Hf.xelem(46) = 0;
          Hf.xelem(47) = 0;
          Hf.xelem(48) = 0;
          Hf.xelem(49) = (r*(r+1)*(1-s2))/2.0E+0;
          Hf.xelem(50) = 0;
          Hf.xelem(51) = 0;
          Hf.xelem(52) = 0;
          Hf.xelem(53) = (r*(r+1)*(1-s2))/2.0E+0;
          Hf.xelem(54) = ((1-r2)*s*(s+1))/2.0E+0;
          Hf.xelem(55) = 0;
          Hf.xelem(56) = 0;
          Hf.xelem(57) = 0;
          Hf.xelem(58) = ((1-r2)*s*(s+1))/2.0E+0;
          Hf.xelem(59) = 0;
          Hf.xelem(60) = 0;
          Hf.xelem(61) = 0;
          Hf.xelem(62) = ((1-r2)*s*(s+1))/2.0E+0;
          Hf.xelem(63) = ((r-1)*r*(1-s2))/2.0E+0;
          Hf.xelem(64) = 0;
          Hf.xelem(65) = 0;
          Hf.xelem(66) = 0;
          Hf.xelem(67) = ((r-1)*r*(1-s2))/2.0E+0;
          Hf.xelem(68) = 0;
          Hf.xelem(69) = 0;
          Hf.xelem(70) = 0;
          Hf.xelem(71) = ((r-1)*r*(1-s2))/2.0E+0;
          Hf.xelem(72) = (1-r2)*(1-s2);
          Hf.xelem(73) = 0;
          Hf.xelem(74) = 0;
          Hf.xelem(75) = 0;
          Hf.xelem(76) = (1-r2)*(1-s2);
          Hf.xelem(77) = 0;
          Hf.xelem(78) = 0;
          Hf.xelem(79) = 0;
          Hf.xelem(80) = (1-r2)*(1-s2);
     }

     static void VectorInterpMatrixDerR(const ColumnVector& rv, Matrix& dHf_dr) {
          FEM_ASSERT(rv.rows() == 2);

          const double r = rv.xelem(0);
          const double s = rv.xelem(1);
          const double s2 = s * s;

          dHf_dr.xelem(0) = (r*(s-1)*s)/4.0E+0+((r-1)*(s-1)*s)/4.0E+0;
          dHf_dr.xelem(1) = 0;
          dHf_dr.xelem(2) = 0;
          dHf_dr.xelem(3) = 0;
          dHf_dr.xelem(4) = (r*(s-1)*s)/4.0E+0+((r-1)*(s-1)*s)/4.0E+0;
          dHf_dr.xelem(5) = 0;
          dHf_dr.xelem(6) = 0;
          dHf_dr.xelem(7) = 0;
          dHf_dr.xelem(8) = (r*(s-1)*s)/4.0E+0+((r-1)*(s-1)*s)/4.0E+0;
          dHf_dr.xelem(9) = ((r+1)*(s-1)*s)/4.0E+0+(r*(s-1)*s)/4.0E+0;
          dHf_dr.xelem(10) = 0;
          dHf_dr.xelem(11) = 0;
          dHf_dr.xelem(12) = 0;
          dHf_dr.xelem(13) = ((r+1)*(s-1)*s)/4.0E+0+(r*(s-1)*s)/4.0E+0;
          dHf_dr.xelem(14) = 0;
          dHf_dr.xelem(15) = 0;
          dHf_dr.xelem(16) = 0;
          dHf_dr.xelem(17) = ((r+1)*(s-1)*s)/4.0E+0+(r*(s-1)*s)/4.0E+0;
          dHf_dr.xelem(18) = ((r+1)*s*(s+1))/4.0E+0+(r*s*(s+1))/4.0E+0;
          dHf_dr.xelem(19) = 0;
          dHf_dr.xelem(20) = 0;
          dHf_dr.xelem(21) = 0;
          dHf_dr.xelem(22) = ((r+1)*s*(s+1))/4.0E+0+(r*s*(s+1))/4.0E+0;
          dHf_dr.xelem(23) = 0;
          dHf_dr.xelem(24) = 0;
          dHf_dr.xelem(25) = 0;
          dHf_dr.xelem(26) = ((r+1)*s*(s+1))/4.0E+0+(r*s*(s+1))/4.0E+0;
          dHf_dr.xelem(27) = (r*s*(s+1))/4.0E+0+((r-1)*s*(s+1))/4.0E+0;
          dHf_dr.xelem(28) = 0;
          dHf_dr.xelem(29) = 0;
          dHf_dr.xelem(30) = 0;
          dHf_dr.xelem(31) = (r*s*(s+1))/4.0E+0+((r-1)*s*(s+1))/4.0E+0;
          dHf_dr.xelem(32) = 0;
          dHf_dr.xelem(33) = 0;
          dHf_dr.xelem(34) = 0;
          dHf_dr.xelem(35) = (r*s*(s+1))/4.0E+0+((r-1)*s*(s+1))/4.0E+0;
          dHf_dr.xelem(36) = -r*(s-1)*s;
          dHf_dr.xelem(37) = 0;
          dHf_dr.xelem(38) = 0;
          dHf_dr.xelem(39) = 0;
          dHf_dr.xelem(40) = -r*(s-1)*s;
          dHf_dr.xelem(41) = 0;
          dHf_dr.xelem(42) = 0;
          dHf_dr.xelem(43) = 0;
          dHf_dr.xelem(44) = -r*(s-1)*s;
          dHf_dr.xelem(45) = ((r+1)*(1-s2))/2.0E+0+(r*(1-s2))/2.0E+0;
          dHf_dr.xelem(46) = 0;
          dHf_dr.xelem(47) = 0;
          dHf_dr.xelem(48) = 0;
          dHf_dr.xelem(49) = ((r+1)*(1-s2))/2.0E+0+(r*(1-s2))/2.0E+0;
          dHf_dr.xelem(50) = 0;
          dHf_dr.xelem(51) = 0;
          dHf_dr.xelem(52) = 0;
          dHf_dr.xelem(53) = ((r+1)*(1-s2))/2.0E+0+(r*(1-s2))/2.0E+0;
          dHf_dr.xelem(54) = -r*s*(s+1);
          dHf_dr.xelem(55) = 0;
          dHf_dr.xelem(56) = 0;
          dHf_dr.xelem(57) = 0;
          dHf_dr.xelem(58) = -r*s*(s+1);
          dHf_dr.xelem(59) = 0;
          dHf_dr.xelem(60) = 0;
          dHf_dr.xelem(61) = 0;
          dHf_dr.xelem(62) = -r*s*(s+1);
          dHf_dr.xelem(63) = (r*(1-s2))/2.0E+0+((r-1)*(1-s2))/2.0E+0;
          dHf_dr.xelem(64) = 0;
          dHf_dr.xelem(65) = 0;
          dHf_dr.xelem(66) = 0;
          dHf_dr.xelem(67) = (r*(1-s2))/2.0E+0+((r-1)*(1-s2))/2.0E+0;
          dHf_dr.xelem(68) = 0;
          dHf_dr.xelem(69) = 0;
          dHf_dr.xelem(70) = 0;
          dHf_dr.xelem(71) = (r*(1-s2))/2.0E+0+((r-1)*(1-s2))/2.0E+0;
          dHf_dr.xelem(72) = -2*r*(1-s2);
          dHf_dr.xelem(73) = 0;
          dHf_dr.xelem(74) = 0;
          dHf_dr.xelem(75) = 0;
          dHf_dr.xelem(76) = -2*r*(1-s2);
          dHf_dr.xelem(77) = 0;
          dHf_dr.xelem(78) = 0;
          dHf_dr.xelem(79) = 0;
          dHf_dr.xelem(80) = -2*r*(1-s2);
     }

     static void VectorInterpMatrixDerS(const ColumnVector& rv, Matrix& dHf_ds) {
          FEM_ASSERT(rv.rows() == 2);

          const double r = rv.xelem(0);
          const double s = rv.xelem(1);
          const double r2 = r * r;

          dHf_ds.xelem(0) = ((r-1)*r*s)/4.0E+0+((r-1)*r*(s-1))/4.0E+0;
          dHf_ds.xelem(1) = 0;
          dHf_ds.xelem(2) = 0;
          dHf_ds.xelem(3) = 0;
          dHf_ds.xelem(4) = ((r-1)*r*s)/4.0E+0+((r-1)*r*(s-1))/4.0E+0;
          dHf_ds.xelem(5) = 0;
          dHf_ds.xelem(6) = 0;
          dHf_ds.xelem(7) = 0;
          dHf_ds.xelem(8) = ((r-1)*r*s)/4.0E+0+((r-1)*r*(s-1))/4.0E+0;
          dHf_ds.xelem(9) = (r*(r+1)*s)/4.0E+0+(r*(r+1)*(s-1))/4.0E+0;
          dHf_ds.xelem(10) = 0;
          dHf_ds.xelem(11) = 0;
          dHf_ds.xelem(12) = 0;
          dHf_ds.xelem(13) = (r*(r+1)*s)/4.0E+0+(r*(r+1)*(s-1))/4.0E+0;
          dHf_ds.xelem(14) = 0;
          dHf_ds.xelem(15) = 0;
          dHf_ds.xelem(16) = 0;
          dHf_ds.xelem(17) = (r*(r+1)*s)/4.0E+0+(r*(r+1)*(s-1))/4.0E+0;
          dHf_ds.xelem(18) = (r*(r+1)*(s+1))/4.0E+0+(r*(r+1)*s)/4.0E+0;
          dHf_ds.xelem(19) = 0;
          dHf_ds.xelem(20) = 0;
          dHf_ds.xelem(21) = 0;
          dHf_ds.xelem(22) = (r*(r+1)*(s+1))/4.0E+0+(r*(r+1)*s)/4.0E+0;
          dHf_ds.xelem(23) = 0;
          dHf_ds.xelem(24) = 0;
          dHf_ds.xelem(25) = 0;
          dHf_ds.xelem(26) = (r*(r+1)*(s+1))/4.0E+0+(r*(r+1)*s)/4.0E+0;
          dHf_ds.xelem(27) = ((r-1)*r*(s+1))/4.0E+0+((r-1)*r*s)/4.0E+0;
          dHf_ds.xelem(28) = 0;
          dHf_ds.xelem(29) = 0;
          dHf_ds.xelem(30) = 0;
          dHf_ds.xelem(31) = ((r-1)*r*(s+1))/4.0E+0+((r-1)*r*s)/4.0E+0;
          dHf_ds.xelem(32) = 0;
          dHf_ds.xelem(33) = 0;
          dHf_ds.xelem(34) = 0;
          dHf_ds.xelem(35) = ((r-1)*r*(s+1))/4.0E+0+((r-1)*r*s)/4.0E+0;
          dHf_ds.xelem(36) = ((1-r2)*s)/2.0E+0+((1-r2)*(s-1))/2.0E+0;
          dHf_ds.xelem(37) = 0;
          dHf_ds.xelem(38) = 0;
          dHf_ds.xelem(39) = 0;
          dHf_ds.xelem(40) = ((1-r2)*s)/2.0E+0+((1-r2)*(s-1))/2.0E+0;
          dHf_ds.xelem(41) = 0;
          dHf_ds.xelem(42) = 0;
          dHf_ds.xelem(43) = 0;
          dHf_ds.xelem(44) = ((1-r2)*s)/2.0E+0+((1-r2)*(s-1))/2.0E+0;
          dHf_ds.xelem(45) = -r*(r+1)*s;
          dHf_ds.xelem(46) = 0;
          dHf_ds.xelem(47) = 0;
          dHf_ds.xelem(48) = 0;
          dHf_ds.xelem(49) = -r*(r+1)*s;
          dHf_ds.xelem(50) = 0;
          dHf_ds.xelem(51) = 0;
          dHf_ds.xelem(52) = 0;
          dHf_ds.xelem(53) = -r*(r+1)*s;
          dHf_ds.xelem(54) = ((1-r2)*(s+1))/2.0E+0+((1-r2)*s)/2.0E+0;
          dHf_ds.xelem(55) = 0;
          dHf_ds.xelem(56) = 0;
          dHf_ds.xelem(57) = 0;
          dHf_ds.xelem(58) = ((1-r2)*(s+1))/2.0E+0+((1-r2)*s)/2.0E+0;
          dHf_ds.xelem(59) = 0;
          dHf_ds.xelem(60) = 0;
          dHf_ds.xelem(61) = 0;
          dHf_ds.xelem(62) = ((1-r2)*(s+1))/2.0E+0+((1-r2)*s)/2.0E+0;
          dHf_ds.xelem(63) = -(r-1)*r*s;
          dHf_ds.xelem(64) = 0;
          dHf_ds.xelem(65) = 0;
          dHf_ds.xelem(66) = 0;
          dHf_ds.xelem(67) = -(r-1)*r*s;
          dHf_ds.xelem(68) = 0;
          dHf_ds.xelem(69) = 0;
          dHf_ds.xelem(70) = 0;
          dHf_ds.xelem(71) = -(r-1)*r*s;
          dHf_ds.xelem(72) = -2*(1-r2)*s;
          dHf_ds.xelem(73) = 0;
          dHf_ds.xelem(74) = 0;
          dHf_ds.xelem(75) = 0;
          dHf_ds.xelem(76) = -2*(1-r2)*s;
          dHf_ds.xelem(77) = 0;
          dHf_ds.xelem(78) = 0;
          dHf_ds.xelem(79) = 0;
          dHf_ds.xelem(80) = -2*(1-r2)*s;
     }

     static void AllocIntegrationRule(Element::FemMatrixType eMatType) {
          static constexpr octave_idx_type N = 3;
          static constexpr double r[2][N] = {{0.774596669241483, 0., -0.774596669241483}, {1., 0., -1.}};
          static constexpr double alpha[2][N] = {{0.555555555555556, 0.888888888888889, 0.555555555555556}, {2./3., 2./3., 2./3.}};

          const octave_idx_type iIntegRule = SelectIntegrationRule(eMatType);

          if (!rgIntegRule[iIntegRule].iGetNumEvalPoints()) {
               rgIntegRule[iIntegRule].SetNumEvalPoints(N * N, 2);

               octave_idx_type l = 0;

               for (octave_idx_type i = 0; i < N; ++i) {
                    for (octave_idx_type j = 0; j < N; ++j) {
                         rgIntegRule[iIntegRule].SetPosition(l, 0, r[iIntegRule][i]);
                         rgIntegRule[iIntegRule].SetPosition(l, 1, r[iIntegRule][j]);
                         rgIntegRule[iIntegRule].SetWeight(l, alpha[iIntegRule][i] * alpha[iIntegRule][j]);
                         ++l;
                    }
               }
          }
     }

     static const IntegrationRule& GetIntegrationRule(Element::FemMatrixType eMatType) {
          const octave_idx_type iIntegRule = SelectIntegrationRule(eMatType);

          FEM_ASSERT(iIntegRule >= 0);
          FEM_ASSERT(static_cast<size_t>(iIntegRule) < rgIntegRule.size());
          FEM_ASSERT(rgIntegRule[iIntegRule].iGetNumEvalPoints() > 0);

          return rgIntegRule[iIntegRule];
     }

private:
    static octave_idx_type SelectIntegrationRule(Element::FemMatrixType eMatType) {
          switch (eMatType) {
          case Element::VEC_LOAD_LUMPED:
               return 1;
          default:
               return 0;
          }
     }

     static array<IntegrationRule, 2> rgIntegRule;
};

array<IntegrationRule, 2> ShapeQuad9::rgIntegRule;

class ShapeTria6H {
public:
     static constexpr octave_idx_type iGetNumNodes() {
          return 6;
     }

     static constexpr octave_idx_type iGetNumDirections() {
          return 2;
     }

     static constexpr octave_idx_type iGetNumDofNode() {
          return 3;
     }

     static constexpr octave_idx_type iGetNumEqualityConstr() {
          return 0;
     }

     static constexpr octave_idx_type iGetNumInequalityConstr() {
          return 1;
     }

     static constexpr double EqualityConstr(const ColumnVector&) {
          return 0;
     }

     static double InequalityConstr(const ColumnVector& rv) {
       return rv.xelem(0) + rv.xelem(1) - 1.;
     }

     static void GetElemLimits(ColumnVector& rmin, ColumnVector& rmax) {
          FEM_ASSERT(rmin.rows() == iGetNumDirections());
          FEM_ASSERT(rmax.rows() == rmin.rows());

          for (octave_idx_type i = 0; i < rmin.rows(); ++i) {
               rmin.xelem(i) = 0.;
          }

          for (octave_idx_type i = 0; i < rmax.rows(); ++i) {
               rmax.xelem(i) = 1.;
          }
     }

     static void ScalarInterpMatrix(const ColumnVector& rv, Matrix& HA, octave_idx_type irow) {
          FEM_ASSERT(rv.rows() == 2);
          FEM_ASSERT(HA.rows() > irow);
          FEM_ASSERT(HA.columns() == 6);

          const double zeta = rv.xelem(0);
          const double eta = rv.xelem(1);
          const octave_idx_type nrows = HA.rows();

          HA.xelem(irow + nrows * 0) = (1-2*((-zeta)-eta+1))*(zeta+eta-1);
          HA.xelem(irow + nrows * 1) = -(1-2*zeta)*zeta;
          HA.xelem(irow + nrows * 2) = -(1-2*eta)*eta;
          HA.xelem(irow + nrows * 3) = 4*((-zeta)-eta+1)*zeta;
          HA.xelem(irow + nrows * 4) = 4*eta*zeta;
          HA.xelem(irow + nrows * 5) = 4*eta*((-zeta)-eta+1);
     }

     static void VectorInterpMatrix(const ColumnVector& rv, Matrix& Hf) {
          FEM_ASSERT(rv.rows() == 2);

          const double zeta = rv.xelem(0);
          const double eta = rv.xelem(1);

          Hf.xelem(0) = (1-2*((-zeta)-eta+1))*(zeta+eta-1);
          Hf.xelem(1) = 0;
          Hf.xelem(2) = 0;
          Hf.xelem(3) = 0;
          Hf.xelem(4) = (1-2*((-zeta)-eta+1))*(zeta+eta-1);
          Hf.xelem(5) = 0;
          Hf.xelem(6) = 0;
          Hf.xelem(7) = 0;
          Hf.xelem(8) = (1-2*((-zeta)-eta+1))*(zeta+eta-1);
          Hf.xelem(9) = -(1-2*zeta)*zeta;
          Hf.xelem(10) = 0;
          Hf.xelem(11) = 0;
          Hf.xelem(12) = 0;
          Hf.xelem(13) = -(1-2*zeta)*zeta;
          Hf.xelem(14) = 0;
          Hf.xelem(15) = 0;
          Hf.xelem(16) = 0;
          Hf.xelem(17) = -(1-2*zeta)*zeta;
          Hf.xelem(18) = -(1-2*eta)*eta;
          Hf.xelem(19) = 0;
          Hf.xelem(20) = 0;
          Hf.xelem(21) = 0;
          Hf.xelem(22) = -(1-2*eta)*eta;
          Hf.xelem(23) = 0;
          Hf.xelem(24) = 0;
          Hf.xelem(25) = 0;
          Hf.xelem(26) = -(1-2*eta)*eta;
          Hf.xelem(27) = 4*((-zeta)-eta+1)*zeta;
          Hf.xelem(28) = 0;
          Hf.xelem(29) = 0;
          Hf.xelem(30) = 0;
          Hf.xelem(31) = 4*((-zeta)-eta+1)*zeta;
          Hf.xelem(32) = 0;
          Hf.xelem(33) = 0;
          Hf.xelem(34) = 0;
          Hf.xelem(35) = 4*((-zeta)-eta+1)*zeta;
          Hf.xelem(36) = 4*eta*zeta;
          Hf.xelem(37) = 0;
          Hf.xelem(38) = 0;
          Hf.xelem(39) = 0;
          Hf.xelem(40) = 4*eta*zeta;
          Hf.xelem(41) = 0;
          Hf.xelem(42) = 0;
          Hf.xelem(43) = 0;
          Hf.xelem(44) = 4*eta*zeta;
          Hf.xelem(45) = 4*eta*((-zeta)-eta+1);
          Hf.xelem(46) = 0;
          Hf.xelem(47) = 0;
          Hf.xelem(48) = 0;
          Hf.xelem(49) = 4*eta*((-zeta)-eta+1);
          Hf.xelem(50) = 0;
          Hf.xelem(51) = 0;
          Hf.xelem(52) = 0;
          Hf.xelem(53) = 4*eta*((-zeta)-eta+1);
     }

     static void VectorInterpMatrixDerR(const ColumnVector& rv, Matrix& dHf_dr) {
          FEM_ASSERT(rv.rows() == 2);

          const double zeta = rv.xelem(0);
          const double eta = rv.xelem(1);

          dHf_dr.xelem(0) = 2*(zeta+eta-1)-2*((-zeta)-eta+1)+1;
          dHf_dr.xelem(1) = 0;
          dHf_dr.xelem(2) = 0;
          dHf_dr.xelem(3) = 0;
          dHf_dr.xelem(4) = 2*(zeta+eta-1)-2*((-zeta)-eta+1)+1;
          dHf_dr.xelem(5) = 0;
          dHf_dr.xelem(6) = 0;
          dHf_dr.xelem(7) = 0;
          dHf_dr.xelem(8) = 2*(zeta+eta-1)-2*((-zeta)-eta+1)+1;
          dHf_dr.xelem(9) = 4*zeta-1;
          dHf_dr.xelem(10) = 0;
          dHf_dr.xelem(11) = 0;
          dHf_dr.xelem(12) = 0;
          dHf_dr.xelem(13) = 4*zeta-1;
          dHf_dr.xelem(14) = 0;
          dHf_dr.xelem(15) = 0;
          dHf_dr.xelem(16) = 0;
          dHf_dr.xelem(17) = 4*zeta-1;
          dHf_dr.xelem(18) = 0;
          dHf_dr.xelem(19) = 0;
          dHf_dr.xelem(20) = 0;
          dHf_dr.xelem(21) = 0;
          dHf_dr.xelem(22) = 0;
          dHf_dr.xelem(23) = 0;
          dHf_dr.xelem(24) = 0;
          dHf_dr.xelem(25) = 0;
          dHf_dr.xelem(26) = 0;
          dHf_dr.xelem(27) = 4*((-zeta)-eta+1)-4*zeta;
          dHf_dr.xelem(28) = 0;
          dHf_dr.xelem(29) = 0;
          dHf_dr.xelem(30) = 0;
          dHf_dr.xelem(31) = 4*((-zeta)-eta+1)-4*zeta;
          dHf_dr.xelem(32) = 0;
          dHf_dr.xelem(33) = 0;
          dHf_dr.xelem(34) = 0;
          dHf_dr.xelem(35) = 4*((-zeta)-eta+1)-4*zeta;
          dHf_dr.xelem(36) = 4*eta;
          dHf_dr.xelem(37) = 0;
          dHf_dr.xelem(38) = 0;
          dHf_dr.xelem(39) = 0;
          dHf_dr.xelem(40) = 4*eta;
          dHf_dr.xelem(41) = 0;
          dHf_dr.xelem(42) = 0;
          dHf_dr.xelem(43) = 0;
          dHf_dr.xelem(44) = 4*eta;
          dHf_dr.xelem(45) = -4*eta;
          dHf_dr.xelem(46) = 0;
          dHf_dr.xelem(47) = 0;
          dHf_dr.xelem(48) = 0;
          dHf_dr.xelem(49) = -4*eta;
          dHf_dr.xelem(50) = 0;
          dHf_dr.xelem(51) = 0;
          dHf_dr.xelem(52) = 0;
          dHf_dr.xelem(53) = -4*eta;
     }

     static void VectorInterpMatrixDerS(const ColumnVector& rv, Matrix& dHf_ds) {
          FEM_ASSERT(rv.rows() == 2);

          const double zeta = rv.xelem(0);
          const double eta = rv.xelem(1);

          dHf_ds.xelem(0) = 2*(zeta+eta-1)-2*((-zeta)-eta+1)+1;
          dHf_ds.xelem(1) = 0;
          dHf_ds.xelem(2) = 0;
          dHf_ds.xelem(3) = 0;
          dHf_ds.xelem(4) = 2*(zeta+eta-1)-2*((-zeta)-eta+1)+1;
          dHf_ds.xelem(5) = 0;
          dHf_ds.xelem(6) = 0;
          dHf_ds.xelem(7) = 0;
          dHf_ds.xelem(8) = 2*(zeta+eta-1)-2*((-zeta)-eta+1)+1;
          dHf_ds.xelem(9) = 0;
          dHf_ds.xelem(10) = 0;
          dHf_ds.xelem(11) = 0;
          dHf_ds.xelem(12) = 0;
          dHf_ds.xelem(13) = 0;
          dHf_ds.xelem(14) = 0;
          dHf_ds.xelem(15) = 0;
          dHf_ds.xelem(16) = 0;
          dHf_ds.xelem(17) = 0;
          dHf_ds.xelem(18) = 4*eta-1;
          dHf_ds.xelem(19) = 0;
          dHf_ds.xelem(20) = 0;
          dHf_ds.xelem(21) = 0;
          dHf_ds.xelem(22) = 4*eta-1;
          dHf_ds.xelem(23) = 0;
          dHf_ds.xelem(24) = 0;
          dHf_ds.xelem(25) = 0;
          dHf_ds.xelem(26) = 4*eta-1;
          dHf_ds.xelem(27) = -4*zeta;
          dHf_ds.xelem(28) = 0;
          dHf_ds.xelem(29) = 0;
          dHf_ds.xelem(30) = 0;
          dHf_ds.xelem(31) = -4*zeta;
          dHf_ds.xelem(32) = 0;
          dHf_ds.xelem(33) = 0;
          dHf_ds.xelem(34) = 0;
          dHf_ds.xelem(35) = -4*zeta;
          dHf_ds.xelem(36) = 4*zeta;
          dHf_ds.xelem(37) = 0;
          dHf_ds.xelem(38) = 0;
          dHf_ds.xelem(39) = 0;
          dHf_ds.xelem(40) = 4*zeta;
          dHf_ds.xelem(41) = 0;
          dHf_ds.xelem(42) = 0;
          dHf_ds.xelem(43) = 0;
          dHf_ds.xelem(44) = 4*zeta;
          dHf_ds.xelem(45) = 4*((-zeta)-eta+1)-4*eta;
          dHf_ds.xelem(46) = 0;
          dHf_ds.xelem(47) = 0;
          dHf_ds.xelem(48) = 0;
          dHf_ds.xelem(49) = 4*((-zeta)-eta+1)-4*eta;
          dHf_ds.xelem(50) = 0;
          dHf_ds.xelem(51) = 0;
          dHf_ds.xelem(52) = 0;
          dHf_ds.xelem(53) = 4*((-zeta)-eta+1)-4*eta;
     }

     static void AllocIntegrationRule(Element::FemMatrixType eMatType) {
          const octave_idx_type iIntegRule = SelectIntegrationRule(eMatType);

          if (!rgIntegRule[iIntegRule].iGetNumEvalPoints()) {
               constexpr double A = 0.470142064105115;
               constexpr double B = 0.101286507323456;
               constexpr double P1 = 0.066197076394253;
               constexpr double P2 = 0.062969590272413;
               static constexpr octave_idx_type N[2] = {6, 7};
               static constexpr double zeta[2][7] = {{0., 1., 0., 1./2., 1./2., 0.}, {1./3., A, 1. - 2. * A, A, B, 1. - 2. * B, B}};
               static constexpr double eta[2][7] = {{0., 0., 1., 0., 1./2., 1./2.}, {1./3., A, A, 1. - 2. * A, B, B, 1. - 2. * B}};
               static constexpr double w[2][7] = {{1./12., 1./12., 1./12., 1./12., 1./12., 1./12.}, {9./80., P1, P1, P1, P2, P2, P2}};

               rgIntegRule[iIntegRule].SetNumEvalPoints(N[iIntegRule], 2);

               for (octave_idx_type i = 0; i < N[iIntegRule]; ++i) {
                    rgIntegRule[iIntegRule].SetPosition(i, 0, zeta[iIntegRule][i]);
                    rgIntegRule[iIntegRule].SetPosition(i, 1, eta[iIntegRule][i]);
                    rgIntegRule[iIntegRule].SetWeight(i, w[iIntegRule][i]);
               }
          }
     }

     static const IntegrationRule& GetIntegrationRule(Element::FemMatrixType eMatType) {
          const octave_idx_type iIntegRule = SelectIntegrationRule(eMatType);

          FEM_ASSERT(iIntegRule >= 0);
          FEM_ASSERT(static_cast<size_t>(iIntegRule) < rgIntegRule.size());
          FEM_ASSERT(rgIntegRule[iIntegRule].iGetNumEvalPoints() > 0);

          return rgIntegRule[iIntegRule];
     }

private:
     static octave_idx_type SelectIntegrationRule(Element::FemMatrixType eMatType) {
          switch (eMatType) {
          case Element::VEC_LOAD_LUMPED:
               return 0;
          default:
               return 1;
          }
     }

     static array<IntegrationRule, 2> rgIntegRule;
};

array<IntegrationRule, 2> ShapeTria6H::rgIntegRule;

class ShapeTria10 {
public:
     static constexpr octave_idx_type iGetNumNodes() {
          return 10;
     }

     static constexpr octave_idx_type iGetNumDirections() {
          return 2;
     }

     static constexpr octave_idx_type iGetNumDofNode() {
          return 3;
     }

     static constexpr octave_idx_type iGetNumEqualityConstr() {
          return 0;
     }

     static constexpr octave_idx_type iGetNumInequalityConstr() {
          return 1;
     }

     static constexpr double EqualityConstr(const ColumnVector&) {
          return 0;
     }

     static double InequalityConstr(const ColumnVector& rv) {
          return rv.xelem(0) + rv.xelem(1) - 1.;
     }

     static void GetElemLimits(ColumnVector& rmin, ColumnVector& rmax) {
          FEM_ASSERT(rmin.rows() == iGetNumDirections());
          FEM_ASSERT(rmax.rows() == rmin.rows());

          for (octave_idx_type i = 0; i < rmin.rows(); ++i) {
               rmin.xelem(i) = 0.;
          }

          for (octave_idx_type i = 0; i < rmax.rows(); ++i) {
               rmax.xelem(i) = 1.;
          }
     }

     static void ScalarInterpMatrix(const ColumnVector& rv, Matrix& HA, octave_idx_type irow) {
          FEM_ASSERT(rv.rows() == 2);
          FEM_ASSERT(HA.rows() > irow);
          FEM_ASSERT(HA.columns() == 10);

          const double zeta = rv.xelem(0);
          const double eta = rv.xelem(1);
          const octave_idx_type nrows = HA.rows();

          HA.xelem(irow + nrows * 0) = (zeta*(3*zeta-2)*(3*zeta-1))/2.0E+0;
          HA.xelem(irow + nrows * 1) = (eta*(3*eta-2)*(3*eta-1))/2.0E+0;
          HA.xelem(irow + nrows * 2) = ((3*((-zeta)-eta+1)-2)*(3*((-zeta)-eta+1)-1)*((-zeta)-eta+1))/2.0E+0;
          HA.xelem(irow + nrows * 3) = (9.0E+0*eta*zeta*(3*zeta-1))/2.0E+0;
          HA.xelem(irow + nrows * 4) = (9.0E+0*eta*(3*eta-1)*zeta)/2.0E+0;
          HA.xelem(irow + nrows * 5) = (9.0E+0*eta*(3*eta-1)*((-zeta)-eta+1))/2.0E+0;
          HA.xelem(irow + nrows * 6) = (9.0E+0*eta*(3*((-zeta)-eta+1)-1)*((-zeta)-eta+1))/2.0E+0;
          HA.xelem(irow + nrows * 7) = (9.0E+0*(3*((-zeta)-eta+1)-1)*((-zeta)-eta+1)*zeta)/2.0E+0;
          HA.xelem(irow + nrows * 8) = (9.0E+0*((-zeta)-eta+1)*zeta*(3*zeta-1))/2.0E+0;
          HA.xelem(irow + nrows * 9) = 27*eta*((-zeta)-eta+1)*zeta;
     }

     static void VectorInterpMatrix(const ColumnVector& rv, Matrix& Hf) {
          FEM_ASSERT(rv.rows() == 2);

          const double zeta = rv.xelem(0);
          const double eta = rv.xelem(1);

          Hf.xelem(0) = (zeta*(3*zeta-2)*(3*zeta-1))/2.0E+0;
          Hf.xelem(1) = 0;
          Hf.xelem(2) = 0;
          Hf.xelem(3) = 0;
          Hf.xelem(4) = (zeta*(3*zeta-2)*(3*zeta-1))/2.0E+0;
          Hf.xelem(5) = 0;
          Hf.xelem(6) = 0;
          Hf.xelem(7) = 0;
          Hf.xelem(8) = (zeta*(3*zeta-2)*(3*zeta-1))/2.0E+0;
          Hf.xelem(9) = (eta*(3*eta-2)*(3*eta-1))/2.0E+0;
          Hf.xelem(10) = 0;
          Hf.xelem(11) = 0;
          Hf.xelem(12) = 0;
          Hf.xelem(13) = (eta*(3*eta-2)*(3*eta-1))/2.0E+0;
          Hf.xelem(14) = 0;
          Hf.xelem(15) = 0;
          Hf.xelem(16) = 0;
          Hf.xelem(17) = (eta*(3*eta-2)*(3*eta-1))/2.0E+0;
          Hf.xelem(18) = ((3*((-zeta)-eta+1)-2)*(3*((-zeta)-eta+1)-1)*((-zeta)-eta+1))/2.0E+0;
          Hf.xelem(19) = 0;
          Hf.xelem(20) = 0;
          Hf.xelem(21) = 0;
          Hf.xelem(22) = ((3*((-zeta)-eta+1)-2)*(3*((-zeta)-eta+1)-1)*((-zeta)-eta+1))/2.0E+0;
          Hf.xelem(23) = 0;
          Hf.xelem(24) = 0;
          Hf.xelem(25) = 0;
          Hf.xelem(26) = ((3*((-zeta)-eta+1)-2)*(3*((-zeta)-eta+1)-1)*((-zeta)-eta+1))/2.0E+0;
          Hf.xelem(27) = (9.0E+0*eta*zeta*(3*zeta-1))/2.0E+0;
          Hf.xelem(28) = 0;
          Hf.xelem(29) = 0;
          Hf.xelem(30) = 0;
          Hf.xelem(31) = (9.0E+0*eta*zeta*(3*zeta-1))/2.0E+0;
          Hf.xelem(32) = 0;
          Hf.xelem(33) = 0;
          Hf.xelem(34) = 0;
          Hf.xelem(35) = (9.0E+0*eta*zeta*(3*zeta-1))/2.0E+0;
          Hf.xelem(36) = (9.0E+0*eta*(3*eta-1)*zeta)/2.0E+0;
          Hf.xelem(37) = 0;
          Hf.xelem(38) = 0;
          Hf.xelem(39) = 0;
          Hf.xelem(40) = (9.0E+0*eta*(3*eta-1)*zeta)/2.0E+0;
          Hf.xelem(41) = 0;
          Hf.xelem(42) = 0;
          Hf.xelem(43) = 0;
          Hf.xelem(44) = (9.0E+0*eta*(3*eta-1)*zeta)/2.0E+0;
          Hf.xelem(45) = (9.0E+0*eta*(3*eta-1)*((-zeta)-eta+1))/2.0E+0;
          Hf.xelem(46) = 0;
          Hf.xelem(47) = 0;
          Hf.xelem(48) = 0;
          Hf.xelem(49) = (9.0E+0*eta*(3*eta-1)*((-zeta)-eta+1))/2.0E+0;
          Hf.xelem(50) = 0;
          Hf.xelem(51) = 0;
          Hf.xelem(52) = 0;
          Hf.xelem(53) = (9.0E+0*eta*(3*eta-1)*((-zeta)-eta+1))/2.0E+0;
          Hf.xelem(54) = (9.0E+0*eta*(3*((-zeta)-eta+1)-1)*((-zeta)-eta+1))/2.0E+0;
          Hf.xelem(55) = 0;
          Hf.xelem(56) = 0;
          Hf.xelem(57) = 0;
          Hf.xelem(58) = (9.0E+0*eta*(3*((-zeta)-eta+1)-1)*((-zeta)-eta+1))/2.0E+0;
          Hf.xelem(59) = 0;
          Hf.xelem(60) = 0;
          Hf.xelem(61) = 0;
          Hf.xelem(62) = (9.0E+0*eta*(3*((-zeta)-eta+1)-1)*((-zeta)-eta+1))/2.0E+0;
          Hf.xelem(63) = (9.0E+0*(3*((-zeta)-eta+1)-1)*((-zeta)-eta+1)*zeta)/2.0E+0;
          Hf.xelem(64) = 0;
          Hf.xelem(65) = 0;
          Hf.xelem(66) = 0;
          Hf.xelem(67) = (9.0E+0*(3*((-zeta)-eta+1)-1)*((-zeta)-eta+1)*zeta)/2.0E+0;
          Hf.xelem(68) = 0;
          Hf.xelem(69) = 0;
          Hf.xelem(70) = 0;
          Hf.xelem(71) = (9.0E+0*(3*((-zeta)-eta+1)-1)*((-zeta)-eta+1)*zeta)/2.0E+0;
          Hf.xelem(72) = (9.0E+0*((-zeta)-eta+1)*zeta*(3*zeta-1))/2.0E+0;
          Hf.xelem(73) = 0;
          Hf.xelem(74) = 0;
          Hf.xelem(75) = 0;
          Hf.xelem(76) = (9.0E+0*((-zeta)-eta+1)*zeta*(3*zeta-1))/2.0E+0;
          Hf.xelem(77) = 0;
          Hf.xelem(78) = 0;
          Hf.xelem(79) = 0;
          Hf.xelem(80) = (9.0E+0*((-zeta)-eta+1)*zeta*(3*zeta-1))/2.0E+0;
          Hf.xelem(81) = 27*eta*((-zeta)-eta+1)*zeta;
          Hf.xelem(82) = 0;
          Hf.xelem(83) = 0;
          Hf.xelem(84) = 0;
          Hf.xelem(85) = 27*eta*((-zeta)-eta+1)*zeta;
          Hf.xelem(86) = 0;
          Hf.xelem(87) = 0;
          Hf.xelem(88) = 0;
          Hf.xelem(89) = 27*eta*((-zeta)-eta+1)*zeta          ;
     }

     static void VectorInterpMatrixDerR(const ColumnVector& rv, Matrix& dHf_dr) {
          FEM_ASSERT(rv.rows() == 2);

          const double zeta = rv.xelem(0);
          const double eta = rv.xelem(1);

          dHf_dr.xelem(0) = ((3*zeta-2)*(3*zeta-1))/2.0E+0+(3.0E+0*zeta*(3*zeta-1))/2.0E+0+(3.0E+0*zeta*(3*zeta-2))/2.0E+0;
          dHf_dr.xelem(1) = 0;
          dHf_dr.xelem(2) = 0;
          dHf_dr.xelem(3) = 0;
          dHf_dr.xelem(4) = ((3*zeta-2)*(3*zeta-1))/2.0E+0+(3.0E+0*zeta*(3*zeta-1))/2.0E+0+(3.0E+0*zeta*(3*zeta-2))/2.0E+0;
          dHf_dr.xelem(5) = 0;
          dHf_dr.xelem(6) = 0;
          dHf_dr.xelem(7) = 0;
          dHf_dr.xelem(8) = ((3*zeta-2)*(3*zeta-1))/2.0E+0+(3.0E+0*zeta*(3*zeta-1))/2.0E+0+(3.0E+0*zeta*(3*zeta-2))/2.0E+0;
          dHf_dr.xelem(9) = 0;
          dHf_dr.xelem(10) = 0;
          dHf_dr.xelem(11) = 0;
          dHf_dr.xelem(12) = 0;
          dHf_dr.xelem(13) = 0;
          dHf_dr.xelem(14) = 0;
          dHf_dr.xelem(15) = 0;
          dHf_dr.xelem(16) = 0;
          dHf_dr.xelem(17) = 0;
          dHf_dr.xelem(18) = ((-3.0E+0)*(3*((-zeta)-eta+1)-1)*((-zeta)-eta+1))/2.0E+0+((-3.0E+0)*(3*((-zeta)-eta+1)-2)*((-zeta)-eta+1))/2.0E+0-((3*((-zeta)-eta+1)-2)*(3*((-zeta)-eta+1)-1))/2.0E+0;
          dHf_dr.xelem(19) = 0;
          dHf_dr.xelem(20) = 0;
          dHf_dr.xelem(21) = 0;
          dHf_dr.xelem(22) = ((-3.0E+0)*(3*((-zeta)-eta+1)-1)*((-zeta)-eta+1))/2.0E+0+((-3.0E+0)*(3*((-zeta)-eta+1)-2)*((-zeta)-eta+1))/2.0E+0-((3*((-zeta)-eta+1)-2)*(3*((-zeta)-eta+1)-1))/2.0E+0;
          dHf_dr.xelem(23) = 0;
          dHf_dr.xelem(24) = 0;
          dHf_dr.xelem(25) = 0;
          dHf_dr.xelem(26) = ((-3.0E+0)*(3*((-zeta)-eta+1)-1)*((-zeta)-eta+1))/2.0E+0+((-3.0E+0)*(3*((-zeta)-eta+1)-2)*((-zeta)-eta+1))/2.0E+0-((3*((-zeta)-eta+1)-2)*(3*((-zeta)-eta+1)-1))/2.0E+0;
          dHf_dr.xelem(27) = (9.0E+0*eta*(3*zeta-1))/2.0E+0+(2.7E+1*eta*zeta)/2.0E+0;
          dHf_dr.xelem(28) = 0;
          dHf_dr.xelem(29) = 0;
          dHf_dr.xelem(30) = 0;
          dHf_dr.xelem(31) = (9.0E+0*eta*(3*zeta-1))/2.0E+0+(2.7E+1*eta*zeta)/2.0E+0;
          dHf_dr.xelem(32) = 0;
          dHf_dr.xelem(33) = 0;
          dHf_dr.xelem(34) = 0;
          dHf_dr.xelem(35) = (9.0E+0*eta*(3*zeta-1))/2.0E+0+(2.7E+1*eta*zeta)/2.0E+0;
          dHf_dr.xelem(36) = (9.0E+0*eta*(3*eta-1))/2.0E+0;
          dHf_dr.xelem(37) = 0;
          dHf_dr.xelem(38) = 0;
          dHf_dr.xelem(39) = 0;
          dHf_dr.xelem(40) = (9.0E+0*eta*(3*eta-1))/2.0E+0;
          dHf_dr.xelem(41) = 0;
          dHf_dr.xelem(42) = 0;
          dHf_dr.xelem(43) = 0;
          dHf_dr.xelem(44) = (9.0E+0*eta*(3*eta-1))/2.0E+0;
          dHf_dr.xelem(45) = ((-9.0E+0)*eta*(3*eta-1))/2.0E+0;
          dHf_dr.xelem(46) = 0;
          dHf_dr.xelem(47) = 0;
          dHf_dr.xelem(48) = 0;
          dHf_dr.xelem(49) = ((-9.0E+0)*eta*(3*eta-1))/2.0E+0;
          dHf_dr.xelem(50) = 0;
          dHf_dr.xelem(51) = 0;
          dHf_dr.xelem(52) = 0;
          dHf_dr.xelem(53) = ((-9.0E+0)*eta*(3*eta-1))/2.0E+0;
          dHf_dr.xelem(54) = ((-2.7E+1)*eta*((-zeta)-eta+1))/2.0E+0+((-9.0E+0)*eta*(3*((-zeta)-eta+1)-1))/2.0E+0;
          dHf_dr.xelem(55) = 0;
          dHf_dr.xelem(56) = 0;
          dHf_dr.xelem(57) = 0;
          dHf_dr.xelem(58) = ((-2.7E+1)*eta*((-zeta)-eta+1))/2.0E+0+((-9.0E+0)*eta*(3*((-zeta)-eta+1)-1))/2.0E+0;
          dHf_dr.xelem(59) = 0;
          dHf_dr.xelem(60) = 0;
          dHf_dr.xelem(61) = 0;
          dHf_dr.xelem(62) = ((-2.7E+1)*eta*((-zeta)-eta+1))/2.0E+0+((-9.0E+0)*eta*(3*((-zeta)-eta+1)-1))/2.0E+0;
          dHf_dr.xelem(63) = ((-2.7E+1)*((-zeta)-eta+1)*zeta)/2.0E+0+((-9.0E+0)*(3*((-zeta)-eta+1)-1)*zeta)/2.0E+0+(9.0E+0*(3*((-zeta)-eta+1)-1)*((-zeta)-eta+1))/2.0E+0;
          dHf_dr.xelem(64) = 0;
          dHf_dr.xelem(65) = 0;
          dHf_dr.xelem(66) = 0;
          dHf_dr.xelem(67) = ((-2.7E+1)*((-zeta)-eta+1)*zeta)/2.0E+0+((-9.0E+0)*(3*((-zeta)-eta+1)-1)*zeta)/2.0E+0+(9.0E+0*(3*((-zeta)-eta+1)-1)*((-zeta)-eta+1))/2.0E+0;
          dHf_dr.xelem(68) = 0;
          dHf_dr.xelem(69) = 0;
          dHf_dr.xelem(70) = 0;
          dHf_dr.xelem(71) = ((-2.7E+1)*((-zeta)-eta+1)*zeta)/2.0E+0+((-9.0E+0)*(3*((-zeta)-eta+1)-1)*zeta)/2.0E+0+(9.0E+0*(3*((-zeta)-eta+1)-1)*((-zeta)-eta+1))/2.0E+0;
          dHf_dr.xelem(72) = ((-9.0E+0)*zeta*(3*zeta-1))/2.0E+0+(9.0E+0*((-zeta)-eta+1)*(3*zeta-1))/2.0E+0+(2.7E+1*((-zeta)-eta+1)*zeta)/2.0E+0;
          dHf_dr.xelem(73) = 0;
          dHf_dr.xelem(74) = 0;
          dHf_dr.xelem(75) = 0;
          dHf_dr.xelem(76) = ((-9.0E+0)*zeta*(3*zeta-1))/2.0E+0+(9.0E+0*((-zeta)-eta+1)*(3*zeta-1))/2.0E+0+(2.7E+1*((-zeta)-eta+1)*zeta)/2.0E+0;
          dHf_dr.xelem(77) = 0;
          dHf_dr.xelem(78) = 0;
          dHf_dr.xelem(79) = 0;
          dHf_dr.xelem(80) = ((-9.0E+0)*zeta*(3*zeta-1))/2.0E+0+(9.0E+0*((-zeta)-eta+1)*(3*zeta-1))/2.0E+0+(2.7E+1*((-zeta)-eta+1)*zeta)/2.0E+0;
          dHf_dr.xelem(81) = 27*eta*((-zeta)-eta+1)-27*eta*zeta;
          dHf_dr.xelem(82) = 0;
          dHf_dr.xelem(83) = 0;
          dHf_dr.xelem(84) = 0;
          dHf_dr.xelem(85) = 27*eta*((-zeta)-eta+1)-27*eta*zeta;
          dHf_dr.xelem(86) = 0;
          dHf_dr.xelem(87) = 0;
          dHf_dr.xelem(88) = 0;
          dHf_dr.xelem(89) = 27*eta*((-zeta)-eta+1)-27*eta*zeta;
     }

     static void VectorInterpMatrixDerS(const ColumnVector& rv, Matrix& dHf_ds) {
          FEM_ASSERT(rv.rows() == 2);

          const double zeta = rv.xelem(0);
          const double eta = rv.xelem(1);

          dHf_ds.xelem(0) = 0;
          dHf_ds.xelem(1) = 0;
          dHf_ds.xelem(2) = 0;
          dHf_ds.xelem(3) = 0;
          dHf_ds.xelem(4) = 0;
          dHf_ds.xelem(5) = 0;
          dHf_ds.xelem(6) = 0;
          dHf_ds.xelem(7) = 0;
          dHf_ds.xelem(8) = 0;
          dHf_ds.xelem(9) = ((3*eta-2)*(3*eta-1))/2.0E+0+(3.0E+0*eta*(3*eta-1))/2.0E+0+(3.0E+0*eta*(3*eta-2))/2.0E+0;
          dHf_ds.xelem(10) = 0;
          dHf_ds.xelem(11) = 0;
          dHf_ds.xelem(12) = 0;
          dHf_ds.xelem(13) = ((3*eta-2)*(3*eta-1))/2.0E+0+(3.0E+0*eta*(3*eta-1))/2.0E+0+(3.0E+0*eta*(3*eta-2))/2.0E+0;
          dHf_ds.xelem(14) = 0;
          dHf_ds.xelem(15) = 0;
          dHf_ds.xelem(16) = 0;
          dHf_ds.xelem(17) = ((3*eta-2)*(3*eta-1))/2.0E+0+(3.0E+0*eta*(3*eta-1))/2.0E+0+(3.0E+0*eta*(3*eta-2))/2.0E+0;
          dHf_ds.xelem(18) = ((-3.0E+0)*(3*((-zeta)-eta+1)-1)*((-zeta)-eta+1))/2.0E+0+((-3.0E+0)*(3*((-zeta)-eta+1)-2)*((-zeta)-eta+1))/2.0E+0-((3*((-zeta)-eta+1)-2)*(3*((-zeta)-eta+1)-1))/2.0E+0;
          dHf_ds.xelem(19) = 0;
          dHf_ds.xelem(20) = 0;
          dHf_ds.xelem(21) = 0;
          dHf_ds.xelem(22) = ((-3.0E+0)*(3*((-zeta)-eta+1)-1)*((-zeta)-eta+1))/2.0E+0+((-3.0E+0)*(3*((-zeta)-eta+1)-2)*((-zeta)-eta+1))/2.0E+0-((3*((-zeta)-eta+1)-2)*(3*((-zeta)-eta+1)-1))/2.0E+0;
          dHf_ds.xelem(23) = 0;
          dHf_ds.xelem(24) = 0;
          dHf_ds.xelem(25) = 0;
          dHf_ds.xelem(26) = ((-3.0E+0)*(3*((-zeta)-eta+1)-1)*((-zeta)-eta+1))/2.0E+0+((-3.0E+0)*(3*((-zeta)-eta+1)-2)*((-zeta)-eta+1))/2.0E+0-((3*((-zeta)-eta+1)-2)*(3*((-zeta)-eta+1)-1))/2.0E+0;
          dHf_ds.xelem(27) = (9.0E+0*zeta*(3*zeta-1))/2.0E+0;
          dHf_ds.xelem(28) = 0;
          dHf_ds.xelem(29) = 0;
          dHf_ds.xelem(30) = 0;
          dHf_ds.xelem(31) = (9.0E+0*zeta*(3*zeta-1))/2.0E+0;
          dHf_ds.xelem(32) = 0;
          dHf_ds.xelem(33) = 0;
          dHf_ds.xelem(34) = 0;
          dHf_ds.xelem(35) = (9.0E+0*zeta*(3*zeta-1))/2.0E+0;
          dHf_ds.xelem(36) = (9.0E+0*(3*eta-1)*zeta)/2.0E+0+(2.7E+1*eta*zeta)/2.0E+0;
          dHf_ds.xelem(37) = 0;
          dHf_ds.xelem(38) = 0;
          dHf_ds.xelem(39) = 0;
          dHf_ds.xelem(40) = (9.0E+0*(3*eta-1)*zeta)/2.0E+0+(2.7E+1*eta*zeta)/2.0E+0;
          dHf_ds.xelem(41) = 0;
          dHf_ds.xelem(42) = 0;
          dHf_ds.xelem(43) = 0;
          dHf_ds.xelem(44) = (9.0E+0*(3*eta-1)*zeta)/2.0E+0+(2.7E+1*eta*zeta)/2.0E+0;
          dHf_ds.xelem(45) = (9.0E+0*(3*eta-1)*((-zeta)-eta+1))/2.0E+0+(2.7E+1*eta*((-zeta)-eta+1))/2.0E+0+((-9.0E+0)*eta*(3*eta-1))/2.0E+0;
          dHf_ds.xelem(46) = 0;
          dHf_ds.xelem(47) = 0;
          dHf_ds.xelem(48) = 0;
          dHf_ds.xelem(49) = (9.0E+0*(3*eta-1)*((-zeta)-eta+1))/2.0E+0+(2.7E+1*eta*((-zeta)-eta+1))/2.0E+0+((-9.0E+0)*eta*(3*eta-1))/2.0E+0;
          dHf_ds.xelem(50) = 0;
          dHf_ds.xelem(51) = 0;
          dHf_ds.xelem(52) = 0;
          dHf_ds.xelem(53) = (9.0E+0*(3*eta-1)*((-zeta)-eta+1))/2.0E+0+(2.7E+1*eta*((-zeta)-eta+1))/2.0E+0+((-9.0E+0)*eta*(3*eta-1))/2.0E+0;
          dHf_ds.xelem(54) = (9.0E+0*(3*((-zeta)-eta+1)-1)*((-zeta)-eta+1))/2.0E+0+((-2.7E+1)*eta*((-zeta)-eta+1))/2.0E+0+((-9.0E+0)*eta*(3*((-zeta)-eta+1)-1))/2.0E+0;
          dHf_ds.xelem(55) = 0;
          dHf_ds.xelem(56) = 0;
          dHf_ds.xelem(57) = 0;
          dHf_ds.xelem(58) = (9.0E+0*(3*((-zeta)-eta+1)-1)*((-zeta)-eta+1))/2.0E+0+((-2.7E+1)*eta*((-zeta)-eta+1))/2.0E+0+((-9.0E+0)*eta*(3*((-zeta)-eta+1)-1))/2.0E+0;
          dHf_ds.xelem(59) = 0;
          dHf_ds.xelem(60) = 0;
          dHf_ds.xelem(61) = 0;
          dHf_ds.xelem(62) = (9.0E+0*(3*((-zeta)-eta+1)-1)*((-zeta)-eta+1))/2.0E+0+((-2.7E+1)*eta*((-zeta)-eta+1))/2.0E+0+((-9.0E+0)*eta*(3*((-zeta)-eta+1)-1))/2.0E+0;
          dHf_ds.xelem(63) = ((-2.7E+1)*((-zeta)-eta+1)*zeta)/2.0E+0+((-9.0E+0)*(3*((-zeta)-eta+1)-1)*zeta)/2.0E+0;
          dHf_ds.xelem(64) = 0;
          dHf_ds.xelem(65) = 0;
          dHf_ds.xelem(66) = 0;
          dHf_ds.xelem(67) = ((-2.7E+1)*((-zeta)-eta+1)*zeta)/2.0E+0+((-9.0E+0)*(3*((-zeta)-eta+1)-1)*zeta)/2.0E+0;
          dHf_ds.xelem(68) = 0;
          dHf_ds.xelem(69) = 0;
          dHf_ds.xelem(70) = 0;
          dHf_ds.xelem(71) = ((-2.7E+1)*((-zeta)-eta+1)*zeta)/2.0E+0+((-9.0E+0)*(3*((-zeta)-eta+1)-1)*zeta)/2.0E+0;
          dHf_ds.xelem(72) = ((-9.0E+0)*zeta*(3*zeta-1))/2.0E+0;
          dHf_ds.xelem(73) = 0;
          dHf_ds.xelem(74) = 0;
          dHf_ds.xelem(75) = 0;
          dHf_ds.xelem(76) = ((-9.0E+0)*zeta*(3*zeta-1))/2.0E+0;
          dHf_ds.xelem(77) = 0;
          dHf_ds.xelem(78) = 0;
          dHf_ds.xelem(79) = 0;
          dHf_ds.xelem(80) = ((-9.0E+0)*zeta*(3*zeta-1))/2.0E+0;
          dHf_ds.xelem(81) = 27*((-zeta)-eta+1)*zeta-27*eta*zeta;
          dHf_ds.xelem(82) = 0;
          dHf_ds.xelem(83) = 0;
          dHf_ds.xelem(84) = 0;
          dHf_ds.xelem(85) = 27*((-zeta)-eta+1)*zeta-27*eta*zeta;
          dHf_ds.xelem(86) = 0;
          dHf_ds.xelem(87) = 0;
          dHf_ds.xelem(88) = 0;
          dHf_ds.xelem(89) = 27*((-zeta)-eta+1)*zeta-27*eta*zeta;
     }

     static void AllocIntegrationRule(Element::FemMatrixType eMatType) {
          const octave_idx_type iIntegRule = SelectIntegrationRule(eMatType);

          static constexpr double zeta2[] = {0.3333333333333330, 0.0288447332326850, 0.7810368490299260, 0.1417072194148800, 0.0250035347626860, 0.0095408154002990};
          static constexpr double eta2[] = {0.3333333333333330, 0.4855776333836570, 0.1094815754850370, 0.3079398387641210, 0.2466725606399030, 0.0668032510122000};
          static constexpr double w2[] = {4.5408995191376998e-02, 1.8362978878233498e-02, 2.2660529717764000e-02, 3.6378958422710002e-02, 1.4163621265528500e-02, 4.7108334818665000e-03};

          static constexpr octave_idx_type perm2[3][3] = {{0, 1, 2}, {1, 0, 2}, {1, 2, 0}};
          static constexpr octave_idx_type perm3[6][3] = {{2, 1, 0}, {2, 0, 1}, {1, 2, 0}, {1, 0, 2}, {0, 2, 1}, {0, 1, 2}};

          if (!rgIntegRule[iIntegRule].iGetNumEvalPoints()) {
               switch (iIntegRule) {
               case R2:
                    rgIntegRule[R2].SetNumEvalPoints(25, 2);
                    rgIntegRule[R2].SetPosition(0, 0, zeta2[0]);
                    rgIntegRule[R2].SetPosition(0, 1, eta2[0]);
                    rgIntegRule[R2].SetWeight(0, w2[0]);

                    for (octave_idx_type i = 1; i < 3; ++i) {
                         for (octave_idx_type j = 0; j < 3; ++j) {
                              const octave_idx_type idx = (i - 1) * 3 + j + 1;

                              FEM_ASSERT(idx >= 1);
                              FEM_ASSERT(idx <= 6);

                              const double x[3] = {zeta2[i], eta2[i], 1. - zeta2[i] - eta2[i]};

                              rgIntegRule[R2].SetWeight(idx, w2[i]);

                              for (octave_idx_type k = 0; k < 2; ++k) {
                                   rgIntegRule[R2].SetPosition(idx, k, x[perm2[j][k]]);
                              }
                         }
                    }

                    for (octave_idx_type i = 3; i < 6; ++i) {
                         for (octave_idx_type j = 0; j < 6; ++j) {
                              const octave_idx_type idx = (i - 3) * 6 + j + 7;

                              FEM_ASSERT(idx >= 7);
                              FEM_ASSERT(idx < 25);

                              const double x[3] = {zeta2[i], eta2[i], 1. - zeta2[i] - eta2[i]};

                              rgIntegRule[R2].SetWeight(idx, w2[i]);

                              for (octave_idx_type k = 0; k < 2; ++k) {
                                   rgIntegRule[R2].SetPosition(idx, k, x[perm3[j][k]]);
                              }
                         }
                    }
                    break;

               default:
                    throw std::logic_error("invalid integration rule");
               }
          }
     }

     static const IntegrationRule& GetIntegrationRule(Element::FemMatrixType eMatType) {
          const octave_idx_type iIntegRule = SelectIntegrationRule(eMatType);

          FEM_ASSERT(iIntegRule >= 0);
          FEM_ASSERT(static_cast<size_t>(iIntegRule) < rgIntegRule.size());
          FEM_ASSERT(rgIntegRule[iIntegRule].iGetNumEvalPoints() > 0);

          return rgIntegRule[iIntegRule];
     }

private:
     enum IntegRuleType {
          R2 = 0,
          RNUM
     };

     static octave_idx_type SelectIntegrationRule(Element::FemMatrixType eMatType) {
          switch (eMatType) {
          case Element::VEC_LOAD_LUMPED:
               throw std::runtime_error("tria10: VEC_LOAD_LUMPED not implemented, use VEC_LOAD_CONSISTENT instead");
          default:
               return R2;
          }
     }

     static array<IntegrationRule, RNUM> rgIntegRule;
};

array<IntegrationRule, ShapeTria10::RNUM> ShapeTria10::rgIntegRule;

class SurfaceElement: public Element {
protected:
     SurfaceElement(ElementTypes::TypeId eltype, octave_idx_type id, const Matrix& X, const Material* material, const int32NDArray& nodes)
          :Element(eltype, id, X, material, nodes) {

          FEM_ASSERT(X.rows() == 3);
          FEM_ASSERT(X.columns() == nodes.numel());
     }

     SurfaceElement(const SurfaceElement& oElem)
          :Element(oElem.eltype, oElem.id, oElem.X, oElem.material, oElem.nodes) {
     }

public:
     virtual void PostProcElem(FemMatrixType eMatType, PostProcData& oSolution) const override {
          switch (eMatType) {
          case VEC_SURFACE_NORMAL_VECTOR:
               SurfaceNormalVectorElem(oSolution.GetField(PostProcData::VEC_EL_SURFACE_NORMAL_VECTOR_RE, eltype), eMatType);
               break;
          case VEC_SURFACE_AREA:
               SurfaceAreaElem(oSolution.GetField(PostProcData::VEC_EL_SURFACE_AREA_RE, eltype), eMatType);
               break;
          default:
               break;
          }
     }

     void SurfaceNormalVectorElem(NDArray& nel, FemMatrixType eMatType) const {
          FEM_ASSERT(nel.ndims() == 3);
          FEM_ASSERT(nel.dim1() >= id);
          FEM_ASSERT(nel.dim2() == nodes.numel());
          FEM_ASSERT(nel.dim3() == 3);

          const IntegrationRule& oIntegRule = GetIntegrationRule(eMatType);
          const octave_idx_type iNumNodes = nodes.numel();
          const octave_idx_type iNumDof = iNumNodes * 3;
          const octave_idx_type iNumDir = oIntegRule.iGetNumDirections();
          const octave_idx_type iNumEvalPoints = oIntegRule.iGetNumEvalPoints();
          const octave_idx_type iNumElem = nel.dim1();

          ColumnVector rv(iNumDir);

          Matrix HA(iNumEvalPoints, iNumNodes);
          ColumnVector n1(3), n2(3), n3(3);
          Matrix dHf_dr(3, iNumDof), dHf_ds(3, iNumDof);
          Matrix ng(iNumEvalPoints, 3);

          for (octave_idx_type i = 0; i < iNumEvalPoints; ++i) {
               for (octave_idx_type j = 0; j < iNumDir; ++j) {
                    rv.xelem(j) = oIntegRule.dGetPosition(i, j);
               }

               ScalarInterpMatrix(rv, HA, i);

               DisplacementInterpMatrixDerR(rv, dHf_dr);
               DisplacementInterpMatrixDerS(rv, dHf_ds);

               SurfaceTangentVector(dHf_dr, n1);
               SurfaceTangentVector(dHf_ds, n2);

               SurfaceNormalVectorUnit(n1, n2, n3);

               for (octave_idx_type j = 0; j < 3; ++j) {
                    ng.xelem(i + iNumEvalPoints * j) = n3.xelem(j);
               }
          }

          const Matrix nn = HA.solve(ng);

          FEM_ASSERT(nn.rows() == iNumNodes);
          FEM_ASSERT(nn.columns() == 3);

          for (octave_idx_type j = 0; j < 3; ++j) {
               for (octave_idx_type i = 0; i < iNumNodes; ++i) {
                    nel.xelem(id - 1 + iNumElem * (i + iNumNodes * j)) = nn.xelem(i + iNumNodes * j);
               }
          }
     }

     void SurfaceAreaElem(NDArray& Ael, FemMatrixType eMatType) const {
          FEM_ASSERT(Ael.ndims() == 2);
          FEM_ASSERT(Ael.dim1() >= id);
          FEM_ASSERT(Ael.dim2() == nodes.numel());

          const IntegrationRule& oIntegRule = GetIntegrationRule(eMatType);
          const octave_idx_type iNumNodes = nodes.numel();
          const octave_idx_type iNumDof = 3 * iNumNodes;
          const octave_idx_type iNumDir = oIntegRule.iGetNumDirections();
          const octave_idx_type iNumElem = Ael.dim1();

          ColumnVector rv(iNumDir);
          Matrix HA(1, iNumNodes);
          ColumnVector n1(3), n2(3);
          Matrix dHf_dr(3, iNumDof), dHf_ds(3, iNumDof);

          for (octave_idx_type i = 0; i < oIntegRule.iGetNumEvalPoints(); ++i) {
               const double alpha = oIntegRule.dGetWeight(i);

               for (octave_idx_type j = 0; j < iNumDir; ++j) {
                    rv.xelem(j) = oIntegRule.dGetPosition(i, j);
               }

               DisplacementInterpMatrixDerR(rv, dHf_dr);
               DisplacementInterpMatrixDerS(rv, dHf_ds);
               ScalarInterpMatrix(rv, HA);

               SurfaceTangentVector(dHf_dr, n1);
               SurfaceTangentVector(dHf_ds, n2);

               const double detJA = JacobianDet(n1, n2);

               for (octave_idx_type m = 0; m < iNumNodes; ++m) {
                    Ael.xelem(id - 1 + iNumElem * m) += alpha * detJA * HA.xelem(m);
               }
          }
     }

     void SurfaceTangentVector(const Matrix& dHf, ColumnVector& n) const {
          FEM_ASSERT(n.rows() == 3);
          FEM_ASSERT(dHf.rows() == 3);
          FEM_ASSERT(dHf.columns() == nodes.numel() * 3);
          FEM_ASSERT(X.rows() == 3);

          for (octave_idx_type i = 0; i < 3; ++i) {
               double ni = 0.;

               for (octave_idx_type j = 0; j < nodes.numel(); ++j) {
                    for (octave_idx_type k = 0; k < 3; ++k) {
                         ni += dHf.xelem(i + 3 * (j * 3 + k)) * X.xelem(k + 3 * j);
                    }
               }

               n.xelem(i) = ni;
          }
     }

     static void SurfaceNormalVector(const ColumnVector& n1, const ColumnVector& n2, ColumnVector& n3) {
          FEM_ASSERT(n1.rows() == 3);
          FEM_ASSERT(n2.rows() == 3);
          FEM_ASSERT(n3.rows() == 3);

          n3.xelem(0) = n1.xelem(1) * n2.xelem(2) - n1.xelem(2) * n2.xelem(1);
          n3.xelem(1) = n1.xelem(2) * n2.xelem(0) - n1.xelem(0) * n2.xelem(2);
          n3.xelem(2) = n1.xelem(0) * n2.xelem(1) - n1.xelem(1) * n2.xelem(0);
     }

     static double SurfaceNormalVectorUnit(const ColumnVector& n1, const ColumnVector& n2, ColumnVector& n3) {
          SurfaceNormalVector(n1, n2, n3);

          const double detJA = JacobianDet(n3);

          for (octave_idx_type i = 0; i < 3; ++i) {
               n3.xelem(i) /= detJA;
          }

          return detJA;
     }

     static double JacobianDet(const ColumnVector& n_detJA) {
          FEM_ASSERT(n_detJA.rows() == 3);

          double detJA_2 = 0.;

          for (octave_idx_type l = 0; l < 3; ++l) {
               detJA_2 += std::pow(n_detJA.xelem(l), 2);
          }

          if (detJA_2 <= 0.) {
               throw std::runtime_error("surface element: Jacobian of surface element is singular");
          }

          return sqrt(detJA_2);
     }

     static double JacobianDet(const ColumnVector& n1, const ColumnVector& n2) {
          FEM_ASSERT(n1.rows() == 3);
          FEM_ASSERT(n2.rows() == 3);

          double detJA = std::pow(n1.xelem(1) * n2.xelem(2) - n1.xelem(2) * n2.xelem(1), 2)
               + std::pow(n1.xelem(2) * n2.xelem(0) - n1.xelem(0) * n2.xelem(2), 2)
               + std::pow(n1.xelem(0) * n2.xelem(1) - n1.xelem(1) * n2.xelem(0), 2);

          if (detJA <= 0.) {
               throw std::runtime_error("surface element: Jacobian of surface element is singular");
          }

          return sqrt(detJA);
     }

     virtual const IntegrationRule& GetIntegrationRule(FemMatrixType eMatType) const=0;
     virtual void ScalarInterpMatrix(const ColumnVector& r, Matrix& HA, octave_idx_type irow = 0) const=0;
     virtual void DisplacementInterpMatrix(const ColumnVector& r, Matrix& Hf) const=0;
     virtual void DisplacementInterpMatrixDerR(const ColumnVector& r, Matrix& dHf_dr) const=0;
     virtual void DisplacementInterpMatrixDerS(const ColumnVector& r, Matrix& dHf_ds) const=0;
};

class PressureElemBase: public SurfaceElement {
public:
protected:
     PressureElemBase(ElementTypes::TypeId eltype, octave_idx_type id, const Matrix& X, const Material* material, const int32NDArray& nodes)
          :SurfaceElement(eltype, id, X, material, nodes) {

          FEM_ASSERT(X.rows() == 3);
          FEM_ASSERT(X.columns() == nodes.numel());
     }

     void AssembleVectors(Matrix& fA, MeshInfo& info, FemMatrixType eMatType, const Matrix& p) const {
          const IntegrationRule& oIntegRule = GetIntegrationRule(eMatType);
          const octave_idx_type iNumNodes = nodes.numel();
          const octave_idx_type iNumDof = iGetNumDof();
          const octave_idx_type iNumDir = oIntegRule.iGetNumDirections();
          const octave_idx_type iNumLoads = p.rows();

          ColumnVector rv(iNumDir);
          Matrix HA(1, iNumNodes), HA_p(1, iNumLoads);
          ColumnVector n1(3), n2(3), n_detJA(3), HfT_n_dA(iNumDof);
          Matrix Hf(3, iNumDof), dHf_dr(3, iNumDof), dHf_ds(3, iNumDof);

          for (octave_idx_type i = 0; i < oIntegRule.iGetNumEvalPoints(); ++i) {
               const double alpha = oIntegRule.dGetWeight(i);

               for (octave_idx_type j = 0; j < iNumDir; ++j) {
                    rv.xelem(j) = oIntegRule.dGetPosition(i, j);
               }

               DisplacementInterpMatrix(rv, Hf);
               DisplacementInterpMatrixDerR(rv, dHf_dr);
               DisplacementInterpMatrixDerS(rv, dHf_ds);
               ScalarInterpMatrix(rv, HA);

               SurfaceTangentVector(dHf_dr, n1);
               SurfaceTangentVector(dHf_ds, n2);
               SurfaceNormalVector(n1, n2, n_detJA);

               info.Add(MeshInfo::JACOBIAN_DET_A, JacobianDet(n_detJA));

               for (octave_idx_type l = 0; l < iNumDof; ++l) {
                    double HfT_nl_detJA = 0.;

                    for (octave_idx_type m = 0; m < 3; ++m) {
                         HfT_nl_detJA -= Hf.xelem(m + 3 * l) * n_detJA.xelem(m);
                    }

                    HfT_n_dA.xelem(l) = HfT_nl_detJA * alpha;
               }

               for (octave_idx_type l = 0; l < iNumLoads; ++l) {
                    double HA_pl = 0.;

                    for (octave_idx_type m = 0; m < iNumNodes; ++m) {
                         HA_pl += HA.xelem(m) * p.xelem(l + iNumLoads * m);
                    }

                    HA_p.xelem(l) = HA_pl;
               }

               for (octave_idx_type l = 0; l < iNumLoads; ++l) {
                    for (octave_idx_type m = 0; m < iNumDof; ++m) {
                         fA.xelem(m + iNumDof * l) += HfT_n_dA.xelem(m) * HA_p.xelem(l);
                    }
               }
          }
     }

     octave_idx_type iGetNumDof() const {
          return nodes.numel() * 3;
     }
};

class PressureLoad: public PressureElemBase {
public:
     PressureLoad(ElementTypes::TypeId eltype, octave_idx_type id, const Matrix& X, const Material* material, const int32NDArray& nodes, const Matrix& p, octave_idx_type colidx)
          :PressureElemBase(eltype, id, X, material, nodes), p(p), colidx(colidx) {

          FEM_ASSERT(X.rows() == 3);
          FEM_ASSERT(X.columns() == p.columns());
          FEM_ASSERT(X.columns() == nodes.numel());
     }

     virtual void Assemble(MatrixAss& mat, MeshInfo& info, const DofMap& dof, FemMatrixType eMatType) const override {
          switch (eMatType) {
          case VEC_LOAD_CONSISTENT:
          case VEC_LOAD_LUMPED:
               break;
          default:
               return;
          }

          const octave_idx_type iNumDof = iGetNumDof();
          const octave_idx_type iNumLoads = p.rows();
          const octave_idx_type iNumNodes = nodes.numel();

          Array<octave_idx_type> dofidx(dim_vector(iNumDof, 1), 0);

          for (octave_idx_type inode = 0; inode < iNumNodes; ++inode) {
               const octave_idx_type inodeidx = nodes.xelem(inode).value() - 1;

               for (octave_idx_type idof = 0; idof < 3; ++idof) {
                    dofidx.xelem(inode * 3 + idof) = dof.GetNodeDofIndex(inodeidx, DofMap::NDOF_DISPLACEMENT, idof);
               }
          }

          Matrix fA(iNumDof, iNumLoads, 0.);

          AssembleVectors(fA, info, eMatType, p);

          for (octave_idx_type j = 0; j < iNumLoads; ++j) {
               for (octave_idx_type i = 0; i < iNumDof; ++i) {
                    mat.Insert(fA.xelem(i + iNumDof * j), dofidx.xelem(i), colidx + j);
               }
          }
     }

     virtual octave_idx_type iGetWorkSpaceSize(FemMatrixType eMatType) const override {
          switch (eMatType) {
          case VEC_LOAD_CONSISTENT:
          case VEC_LOAD_LUMPED:
               return iGetNumDof() * p.rows();
          default:
               return 0;
          }
     }

private:
     const Matrix p;
     const octave_idx_type colidx;
};

class FluidStructInteract: public PressureElemBase {
public:
     FluidStructInteract(ElementTypes::TypeId eltype, octave_idx_type id, const Matrix& X, const Material* material, const int32NDArray& nodes)
          :PressureElemBase(eltype, id, X, material, nodes) {

          FEM_ASSERT(X.rows() == 3);
          FEM_ASSERT(X.columns() == nodes.numel());
     }

     virtual void Assemble(MatrixAss& mat, MeshInfo& info, const DofMap& dof, FemMatrixType eMatType) const override {
          switch (eMatType) {
          case MAT_DAMPING_FLUID_STRUCT_RE:
               break;
          default:
               return;
          }

          const octave_idx_type iNumDofStruct = iGetNumDof();
          const octave_idx_type iNumNodes = nodes.numel();
          const octave_idx_type iNumDofFluid = iNumNodes;

          Array<octave_idx_type> dofidx_s(dim_vector(iNumDofStruct, 1), 0);
          Array<octave_idx_type> dofidx_f(dim_vector(iNumDofFluid, 1), 0);

          for (octave_idx_type inode = 0; inode < iNumNodes; ++inode) {
               const octave_idx_type inodeidx = nodes(inode).value() - 1;

               for (octave_idx_type idof = 0; idof < 3; ++idof) {
                    dofidx_s.xelem(inode * 3 + idof) = dof.GetNodeDofIndex(inodeidx, DofMap::NDOF_DISPLACEMENT, idof);
               }

               dofidx_f.xelem(inode) = dof.GetNodeDofIndex(inodeidx, DofMap::NDOF_VELOCITY_POT, 0);
          }

          Matrix A(iNumDofStruct, iNumDofFluid, 0.);
          Matrix p(iNumDofFluid, iNumDofFluid, 0.);

          for (octave_idx_type i = 0; i < iNumDofFluid; ++i) {
               p.xelem(i + iNumDofFluid * i) = 1.; // Assume, that the mesh is oriented using the solid element volume
          }

          AssembleVectors(A, info, eMatType, p);

          for (octave_idx_type j = 0; j < iNumDofFluid; ++j) {
               for (octave_idx_type i = 0; i < iNumDofStruct; ++i) {
                    const double Aij = A.xelem(i + iNumDofStruct * j);
                    mat.Insert(Aij, dofidx_s.xelem(i), dofidx_f.xelem(j));
                    mat.Insert(Aij, dofidx_f.xelem(j), dofidx_s.xelem(i));
               }
          }
     }

     virtual octave_idx_type iGetWorkSpaceSize(FemMatrixType eMatType) const override {
          switch (eMatType) {
          case MAT_DAMPING_FLUID_STRUCT_RE:
               return 2 * iGetNumDof() * nodes.numel();
          default:
               return 0;
          }
     }
};

class ScalarFieldBC: public SurfaceElement {
public:
     ScalarFieldBC(ElementTypes::TypeId eltype, octave_idx_type id, const Matrix& X, const Material* material, const int32NDArray& nodes)
          :SurfaceElement(eltype, id, X, material, nodes) {

          FEM_ASSERT(X.rows() == 3);
     }

     ScalarFieldBC(const ScalarFieldBC& oElem)=default;

     void CoefficientMatrix(MatrixAss& mat, MeshInfo& info, const DofMap& dof, const FemMatrixType eMatType, const DofMap::NodalDofType eDofType, const RowVector& h) const {
          const IntegrationRule& oIntegRule = GetIntegrationRule(eMatType);
          const octave_idx_type iNumNodes = nodes.numel();
          const octave_idx_type iNumDof = iGetNumDof();
          const octave_idx_type iNumDir = oIntegRule.iGetNumDirections();

          ColumnVector rv(iNumDir);
          Array<octave_idx_type> dofidx(dim_vector(iNumDof, 1), 0);

          for (octave_idx_type inode = 0; inode < iNumNodes; ++inode) {
               dofidx.xelem(inode) = dof.GetNodeDofIndex(nodes(inode).value() - 1, eDofType, 0);
          }

          Matrix HA(1, iNumDof);
          ColumnVector n1(3), n2(3);
          Matrix dHf_dr(3, 3 * iNumDof), dHf_ds(3, 3 * iNumDof), KA(iNumDof, iNumDof, 0.);

          for (octave_idx_type i = 0; i < oIntegRule.iGetNumEvalPoints(); ++i) {
               const double alpha = oIntegRule.dGetWeight(i);

               for (octave_idx_type j = 0; j < iNumDir; ++j) {
                    rv.xelem(j) = oIntegRule.dGetPosition(i, j);
               }

               DisplacementInterpMatrixDerR(rv, dHf_dr);
               DisplacementInterpMatrixDerS(rv, dHf_ds);
               ScalarInterpMatrix(rv, HA);

               SurfaceTangentVector(dHf_dr, n1);
               SurfaceTangentVector(dHf_ds, n2);

               const double detJA = JacobianDet(n1, n2);

               info.Add(MeshInfo::JACOBIAN_DET_A, detJA);

               double hi = 0.;

               for (octave_idx_type j = 0; j < iNumNodes; ++j) {
                    hi += h.xelem(j) * HA.xelem(j);
               }

               hi *= detJA * alpha;

               for (octave_idx_type l = 0; l < iNumDof; ++l) {
                    for (octave_idx_type m = l; m < iNumDof; ++m) {
                         KA.xelem(m + iNumDof * l) += hi * HA.xelem(m) * HA.xelem(l);
                    }
               }
          }

          for (octave_idx_type i = 1; i < iNumDof; ++i) {
               for (octave_idx_type j = 0; j < i; ++j) {
                    KA.xelem(j + iNumDof * i) = KA.xelem(i + iNumDof * j);
               }
          }

          for (octave_idx_type j = 0; j < iNumDof; ++j) {
               for (octave_idx_type i = 0; i < iNumDof; ++i) {
                    mat.Insert(KA.xelem(i + iNumDof * j), dofidx.xelem(i), dofidx.xelem(j));
               }
          }
     }

     void RightHandSideVector(MatrixAss& mat, MeshInfo& info, const DofMap& dof, const FemMatrixType eMatType, const DofMap::NodalDofType eDofType, const Matrix& Thetae, const RowVector& h) const {
          const IntegrationRule& oIntegRule = GetIntegrationRule(eMatType);
          const octave_idx_type iNumNodes = nodes.numel();
          const octave_idx_type iNumDof = iGetNumDof();
          const octave_idx_type iNumDir = oIntegRule.iGetNumDirections();
          const octave_idx_type iNumLoads = Thetae.rows();

          ColumnVector rv(iNumDir);
          Array<octave_idx_type> dofidx(dim_vector(iNumDof, 1), 0);

          for (octave_idx_type inode = 0; inode < iNumNodes; ++inode) {
               dofidx.xelem(inode) = dof.GetNodeDofIndex(nodes(inode).value() - 1, eDofType, 0);
          }

          Matrix HA(1, iNumDof), HA_Thetae(1, iNumLoads);
          ColumnVector n1(3), n2(3);
          Matrix dHf_dr(3, 3 * iNumDof), dHf_ds(3, 3 * iNumDof), QA(iNumDof, iNumLoads, 0.);

          for (octave_idx_type i = 0; i < oIntegRule.iGetNumEvalPoints(); ++i) {
               const double alpha = oIntegRule.dGetWeight(i);

               for (octave_idx_type j = 0; j < iNumDir; ++j) {
                    rv.xelem(j) = oIntegRule.dGetPosition(i, j);
               }

               DisplacementInterpMatrixDerR(rv, dHf_dr);
               DisplacementInterpMatrixDerS(rv, dHf_ds);
               ScalarInterpMatrix(rv, HA);

               SurfaceTangentVector(dHf_dr, n1);
               SurfaceTangentVector(dHf_ds, n2);

               const double detJA = JacobianDet(n1, n2);

               info.Add(MeshInfo::JACOBIAN_DET_A, detJA);

               double hi = 0.;

               for (octave_idx_type j = 0; j < iNumNodes; ++j) {
                    hi += h.xelem(j) * HA.xelem(j);
               }

               for (octave_idx_type l = 0; l < iNumLoads; ++l) {
                    double HA_Thetael = 0.;

                    for (octave_idx_type m = 0; m < iNumNodes; ++m) {
                         HA_Thetael += HA.xelem(m) * Thetae.xelem(l + iNumLoads * m);
                    }

                    HA_Thetae.xelem(l) = HA_Thetael;
               }

               hi *= detJA * alpha;

               for (octave_idx_type l = 0; l < iNumLoads; ++l) {
                    for (octave_idx_type m = 0; m < iNumDof; ++m) {
                         QA.xelem(m + iNumDof * l) += hi * HA.xelem(m) * HA_Thetae.xelem(l);
                    }
               }
          }

          for (octave_idx_type j = 0; j < iNumLoads; ++j) {
               for (octave_idx_type i = 0; i < iNumDof; ++i) {
                    mat.Insert(QA.xelem(i + iNumDof * j), dofidx.xelem(i), j + 1);
               }
          }
     }

     octave_idx_type iGetNumDof() const {
          return nodes.numel();
     }
};

class ThermalConvectionBC: public ScalarFieldBC {
public:
     ThermalConvectionBC(ElementTypes::TypeId eltype, octave_idx_type id, const Matrix& X, const Material* material, const int32NDArray& nodes, const Matrix& Thetae, const RowVector& h)
          :ScalarFieldBC(eltype, id, X, material, nodes), Thetae(Thetae), h(h) {
          FEM_ASSERT(X.columns() == Thetae.columns());
          FEM_ASSERT(X.columns() == nodes.numel());
     }

     ThermalConvectionBC(const ThermalConvectionBC& oElem)=default;

     virtual void Assemble(MatrixAss& mat, MeshInfo& info, const DofMap& dof, const FemMatrixType eMatType) const override final {
          switch (eMatType) {
          case VEC_LOAD_THERMAL:
               RightHandSideVector(mat, info, dof, eMatType, DofMap::NDOF_TEMPERATURE, Thetae, h);
               break;

          case MAT_THERMAL_COND:
               CoefficientMatrix(mat, info, dof, eMatType, DofMap::NDOF_TEMPERATURE, h);
               break;

          default:
               ;
          }
     }

     virtual octave_idx_type iGetWorkSpaceSize(FemMatrixType eMatType) const override final {
          switch (eMatType) {
          case MAT_THERMAL_COND:
               return iGetNumDof() * iGetNumDof();

          case VEC_LOAD_THERMAL:
               return iGetNumDof() * Thetae.rows();

          default:
               return 0;
          }
     }

private:
     const Matrix Thetae;
     const RowVector h;
};

class ParticleVelocityBC: public ScalarFieldBC {
public:
     ParticleVelocityBC(ElementTypes::TypeId eltype, octave_idx_type id, const Matrix& X, const Material* material, const int32NDArray& nodes, const Matrix& vn, const RowVector& coef)
          :ScalarFieldBC(eltype, id, X, material, nodes), vn(vn), coef(coef) {
     }

     ParticleVelocityBC(const ParticleVelocityBC& oElem)=default;

     virtual void Assemble(MatrixAss& mat, MeshInfo& info, const DofMap& dof, const FemMatrixType eMatType) const override final {
          switch (eMatType) {
          case VEC_LOAD_ACOUSTICS:
          case VEC_LOAD_FLUID_STRUCT:
               RightHandSideVector(mat, info, dof, eMatType, DofMap::NDOF_VELOCITY_POT, vn, coef);
               break;

          default:
               ;
          }
     }

     virtual octave_idx_type iGetWorkSpaceSize(FemMatrixType eMatType) const override final {
          switch (eMatType) {
          case VEC_LOAD_ACOUSTICS:
          case VEC_LOAD_FLUID_STRUCT:
               return iGetNumDof() * vn.rows();

          default:
               return 0;
          }
     }

private:
     const Matrix vn;
     const RowVector coef;
};

class AcousticImpedanceBC: public ScalarFieldBC {
public:
     AcousticImpedanceBC(ElementTypes::TypeId eltype, octave_idx_type id, const Matrix& X, const Material* material, const int32NDArray& nodes, const RowVector& rec_z_re, const RowVector& rec_z_im)
          :ScalarFieldBC(eltype, id, X, material, nodes), rec_z_re(rec_z_re), rec_z_im(rec_z_im) {
     }

     AcousticImpedanceBC(const AcousticImpedanceBC& oElem)=default;

     virtual void Assemble(MatrixAss& mat, MeshInfo& info, const DofMap& dof, const FemMatrixType eMatType) const override final {
          switch (eMatType) {
          case MAT_DAMPING_ACOUSTICS_RE:
          case MAT_DAMPING_FLUID_STRUCT_RE:
               CoefficientMatrix(mat, info, dof, eMatType, DofMap::NDOF_VELOCITY_POT, rec_z_re);
               break;
          case MAT_DAMPING_ACOUSTICS_IM:
          case MAT_DAMPING_FLUID_STRUCT_IM:
               CoefficientMatrix(mat, info, dof, eMatType, DofMap::NDOF_VELOCITY_POT, rec_z_im);
               break;

          default:
               ;
          }
     }

     virtual octave_idx_type iGetWorkSpaceSize(FemMatrixType eMatType) const override final {
          switch (eMatType) {
          case MAT_DAMPING_ACOUSTICS_RE:
          case MAT_DAMPING_ACOUSTICS_IM:
          case MAT_DAMPING_FLUID_STRUCT_RE:
          case MAT_DAMPING_FLUID_STRUCT_IM:
               return iGetNumDof() * iGetNumDof();

          default:
               return 0;
          }
     }

private:
     const RowVector rec_z_re, rec_z_im;
};

class AcousticBoundary: public SurfaceElement {
public:
     AcousticBoundary(ElementTypes::TypeId eltype, octave_idx_type id, const Matrix& X, const Material* material, const int32NDArray& nodes)
          :SurfaceElement(eltype, id, X, material, nodes) {

          FEM_ASSERT(X.rows() == 3);
     }

     virtual void Assemble(MatrixAss& mat, MeshInfo& info, const DofMap& dof, const FemMatrixType eMatType) const override final {
          // unused
     }

     virtual void PostProcElem(FemMatrixType eMatType, PostProcData& oSolution) const final override {
          switch (eMatType) {
          case VEC_PARTICLE_VELOCITY:
               ParticleVelocityNormal<double>(oSolution.GetField(PostProcData::SCA_EL_ACOUSTIC_PART_VEL_NORM_RE, eltype),
                                              eMatType,
                                              oSolution.GetField(PostProcData::VEC_NO_ACOUSTIC_PART_VEL_RE, eltype));
               break;
          case VEC_PARTICLE_VELOCITY_C:
               ParticleVelocityNormal<std::complex<double> >(oSolution.GetField(PostProcData::SCA_EL_ACOUSTIC_PART_VEL_NORM_C, eltype),
                                                             eMatType,
                                                             oSolution.GetField(PostProcData::VEC_NO_ACOUSTIC_PART_VEL_C, eltype));
               break;
          case SCA_ACOUSTIC_INTENSITY:
               AcousticIntensity<double>(oSolution.GetField(PostProcData::SCA_EL_ACOUSTIC_INTENSITY_RE, eltype),
                                         oSolution.GetField(PostProcData::SCA_EL_ACOUSTIC_SOUND_POWER_RE, eltype),
                                         eMatType,
                                         oSolution.GetField(PostProcData::VEC_NO_ACOUSTIC_PART_VEL_RE, eltype),
                                         oSolution.GetField(PostProcData::SCA_NO_ACOUSTIC_PART_VEL_POT_P_RE, eltype));
               break;
          case SCA_ACOUSTIC_INTENSITY_C:
               AcousticIntensity<std::complex<double> >(oSolution.GetField(PostProcData::SCA_EL_ACOUSTIC_INTENSITY_RE, eltype),
                                                        oSolution.GetField(PostProcData::SCA_EL_ACOUSTIC_SOUND_POWER_RE, eltype),
                                                        eMatType,
                                                        oSolution.GetField(PostProcData::VEC_NO_ACOUSTIC_PART_VEL_C, eltype),
                                                        oSolution.GetField(PostProcData::SCA_NO_ACOUSTIC_PART_VEL_POT_P_C, eltype));
               break;
          default:
               break;
          }
     }

     template <typename T>
     void ParticleVelocityNormal(typename PostProcTypeTraits<T>::NDArrayType& vn,
                                 FemMatrixType eMatType,
                                 const typename PostProcTypeTraits<T>::NDArrayType& v) const {
          typedef typename PostProcTypeTraits<T>::MatrixType TMatrix;
          typedef typename PostProcTypeTraits<T>::ColumnVectorType TColumnVector;

          const octave_idx_type iNumNodes = nodes.numel();
          const octave_idx_type iNumLoads = v.ndims() >= 3 ? v.dim3() : 1;
          const octave_idx_type iNumNodesMesh = v.dim1();
          const octave_idx_type iNumElem = vn.dim1();

          FEM_ASSERT(v.ndims() >= 2);
          FEM_ASSERT(v.dim2() == 3);
          FEM_ASSERT(vn.ndims() >= 2);
          FEM_ASSERT(vn.dim2() == iNumNodes);
          FEM_ASSERT(vn.ndims() >= 3 ? vn.dim3() == iNumLoads : iNumLoads == 1);

          TMatrix ve(iNumNodes, 3 * iNumLoads);

          for (octave_idx_type k = 0; k < iNumLoads; ++k) {
               for (octave_idx_type j = 0; j < 3; ++j) {
                    for (octave_idx_type i = 0; i < iNumNodes; ++i) {
                         const octave_idx_type inode = nodes.xelem(i).value() - 1;
                         ve.xelem(i + iNumNodes * (k * 3 + j)) = v.xelem(inode + iNumNodesMesh * (j + 3 * k));
                    }
               }
          }

          const IntegrationRule& oIntegRule = GetIntegrationRule(eMatType);
          const octave_idx_type iNumDir = oIntegRule.iGetNumDirections();
          const octave_idx_type iNumEvalPoints = oIntegRule.iGetNumEvalPoints();

          ColumnVector rv(iNumDir);

          Matrix HA(iNumEvalPoints, iNumNodes);
          ColumnVector n1(3), n2(3), n(3);
          TColumnVector vik(3);
          Matrix dHf_dr(3, 3 * iNumNodes), dHf_ds(3, 3 * iNumNodes);
          TMatrix vng(iNumEvalPoints, iNumLoads);

          for (octave_idx_type i = 0; i < iNumEvalPoints; ++i) {
               for (octave_idx_type j = 0; j < iNumDir; ++j) {
                    rv.xelem(j) = oIntegRule.dGetPosition(i, j);
               }

               DisplacementInterpMatrixDerR(rv, dHf_dr);
               DisplacementInterpMatrixDerS(rv, dHf_ds);
               ScalarInterpMatrix(rv, HA, i);

               SurfaceTangentVector(dHf_dr, n1);
               SurfaceTangentVector(dHf_ds, n2);
               SurfaceNormalVectorUnit(n1, n2, n);

               for (octave_idx_type k = 0; k < iNumLoads; ++k) {
                    for (octave_idx_type l = 0; l < 3; ++l) {
                         T vikl{};

                         for (octave_idx_type m = 0; m < iNumNodes; ++m) {
                              vikl += HA.xelem(i + iNumEvalPoints * m) * ve.xelem(m + iNumNodes * (k * 3 + l));
                         }

                         vik.xelem(l) = vikl;
                    }

                    T vnik{};

                    for (octave_idx_type l = 0; l < 3; ++l) {
                         vnik += n.xelem(l) * vik.xelem(l);
                    }

                    vng.xelem(i + iNumEvalPoints * k) = vnik;
               }
          }

          const TMatrix vne = HA.solve(vng);

          FEM_ASSERT(vne.rows() == iNumNodes);
          FEM_ASSERT(vne.columns() == iNumLoads);

          for (octave_idx_type j = 0; j < iNumLoads; ++j) {
               for (octave_idx_type i = 0; i < iNumNodes; ++i) {
                    vn.xelem(id - 1 + iNumElem * (i + iNumNodes * j)) = vne.xelem(i + iNumNodes * j);
               }
          }
     }

     template <typename T>
     void AcousticIntensity(NDArray& I,
                            NDArray& P,
                            FemMatrixType eMatType,
                            const typename PostProcTypeTraits<T>::NDArrayType& v,
                            const typename PostProcTypeTraits<T>::NDArrayType& PhiP) const {
          const octave_idx_type iNumNodes = nodes.numel();
          const octave_idx_type iNumLoads = v.ndims() >= 3 ? v.dim3() : 1;
          const octave_idx_type iNumNodesMesh = v.dim1();
          const octave_idx_type iNumElem = I.dim1();
          typedef typename PostProcTypeTraits<T>::MatrixType TMatrix;
          typedef typename PostProcTypeTraits<T>::ColumnVectorType TColumnVector;

          FEM_ASSERT(v.ndims() >= 2);
          FEM_ASSERT(v.dim2() == 3);
          FEM_ASSERT(I.ndims() >= 2);
          FEM_ASSERT(I.dim2() == iNumNodes);
          FEM_ASSERT(I.ndims() >= 3 ? I.dim3() == iNumLoads : iNumLoads == 1);
          FEM_ASSERT(P.dim2() == iNumLoads);
          FEM_ASSERT(P.dim1() == I.dim1());

          TMatrix ve(iNumNodes, 3 * iNumLoads);

          for (octave_idx_type k = 0; k < iNumLoads; ++k) {
               for (octave_idx_type j = 0; j < 3; ++j) {
                    for (octave_idx_type i = 0; i < iNumNodes; ++i) {
                         const octave_idx_type inode = nodes.xelem(i).value() - 1;
                         ve.xelem(i + iNumNodes * (k * 3 + j)) = v.xelem(inode + iNumNodesMesh * (j + 3 * k));
                    }
               }
          }

          TMatrix pe(iNumNodes, iNumLoads);

          for (octave_idx_type k = 0; k < iNumLoads; ++k) {
               for (octave_idx_type i = 0; i < iNumNodes; ++i) {
                    const octave_idx_type inode = nodes.xelem(i).value() - 1;
                    pe.xelem(i + iNumNodes * k) = -PhiP.xelem(inode + iNumNodesMesh * k);
               }
          }

          const IntegrationRule& oIntegRule = GetIntegrationRule(eMatType);
          const octave_idx_type iNumDir = oIntegRule.iGetNumDirections();
          const octave_idx_type iNumEvalPoints = oIntegRule.iGetNumEvalPoints();

          ColumnVector rv(iNumDir);

          Matrix HA(iNumEvalPoints, iNumNodes);
          ColumnVector n1(3), n2(3), n(3);
          TColumnVector vik(3);
          Matrix dHf_dr(3, 3 * iNumNodes), dHf_ds(3, 3 * iNumNodes);
          Matrix Ig(iNumEvalPoints, iNumLoads);
          RowVector Pe(iNumLoads, 0.);

          for (octave_idx_type i = 0; i < iNumEvalPoints; ++i) {
               const double alphai = oIntegRule.dGetWeight(i);

               for (octave_idx_type j = 0; j < iNumDir; ++j) {
                    rv.xelem(j) = oIntegRule.dGetPosition(i, j);
               }

               DisplacementInterpMatrixDerR(rv, dHf_dr);
               DisplacementInterpMatrixDerS(rv, dHf_ds);
               ScalarInterpMatrix(rv, HA, i);

               SurfaceTangentVector(dHf_dr, n1);
               SurfaceTangentVector(dHf_ds, n2);

               const double detJA = SurfaceNormalVectorUnit(n1, n2, n);

               for (octave_idx_type k = 0; k < iNumLoads; ++k) {
                    for (octave_idx_type l = 0; l < 3; ++l) {
                         T vikl{};

                         for (octave_idx_type m = 0; m < iNumNodes; ++m) {
                              vikl += HA.xelem(i + iNumEvalPoints * m) * ve.xelem(m + iNumNodes * (k * 3 + l));
                         }

                         vik.xelem(l) = vikl;
                    }

                    T vnik{};

                    for (octave_idx_type l = 0; l < 3; ++l) {
                         vnik += n.xelem(l) * vik.xelem(l);
                    }

                    T pik{};

                    for (octave_idx_type m = 0; m < iNumNodes; ++m) {
                         pik += HA.xelem(i + iNumEvalPoints * m) * pe.xelem(m + iNumNodes * k);
                    }

                    const double Iik = PostProcTypeTraits<T>::EffectiveAmplitude(pik, vnik);

                    Ig.xelem(i + iNumEvalPoints * k) = Iik;
                    Pe.xelem(k) += Iik * alphai * detJA;
               }
          }

          const Matrix Ie = HA.solve(Ig);

          FEM_ASSERT(Ie.rows() == iNumNodes);
          FEM_ASSERT(Ie.columns() == iNumLoads);

          for (octave_idx_type j = 0; j < iNumLoads; ++j) {
               for (octave_idx_type i = 0; i < iNumNodes; ++i) {
                    I.xelem(id - 1 + iNumElem * (i + iNumNodes * j)) = Ie.xelem(i + iNumNodes * j);
               }
          }

          for (octave_idx_type j = 0; j < iNumLoads; ++j) {
               P.xelem(id - 1 + iNumElem * j) = Pe.xelem(j);
          }
     }

     virtual octave_idx_type iGetWorkSpaceSize(FemMatrixType eMatType) const override final {
          return 0;
     }
};

class HeatSource: public SurfaceElement {
public:
     HeatSource(ElementTypes::TypeId eltype, octave_idx_type id, const Matrix& X, const Material* material, const int32NDArray& nodes, const Matrix& qe, octave_idx_type colidx)
          :SurfaceElement(eltype, id, X, material, nodes), qe(qe), colidx(colidx) {

          FEM_ASSERT(X.rows() == 3);
          FEM_ASSERT(X.columns() == qe.columns());
          FEM_ASSERT(X.columns() == nodes.numel());
     }

     HeatSource(const HeatSource& oElem)
          :SurfaceElement(oElem.eltype, oElem.id, oElem.X, oElem.material, oElem.nodes), qe(oElem.qe), colidx(oElem.colidx) {
     }

     virtual void Assemble(MatrixAss& mat, MeshInfo& info, const DofMap& dof, const FemMatrixType eMatType) const override {
          switch (eMatType) {
          case VEC_LOAD_THERMAL:
               ThermalLoadVector(mat, info, dof, eMatType);
               break;

          default:
               ;
          }
     }

     void ThermalLoadVector(MatrixAss& mat, MeshInfo& info, const DofMap& dof, const FemMatrixType eMatType) const {
          const IntegrationRule& oIntegRule = GetIntegrationRule(eMatType);
          const octave_idx_type iNumNodes = nodes.numel();
          const octave_idx_type iNumDof = iGetNumDof();
          const octave_idx_type iNumDir = oIntegRule.iGetNumDirections();
          const octave_idx_type iNumLoads = qe.rows();

          ColumnVector rv(iNumDir);
          Array<octave_idx_type> dofidx(dim_vector(iNumDof, 1), 0);

          for (octave_idx_type inode = 0; inode < iNumNodes; ++inode) {
               dofidx.xelem(inode) = dof.GetNodeDofIndex(nodes(inode).value() - 1, DofMap::NDOF_TEMPERATURE, 0);
          }

          Matrix HA(1, iNumDof), HA_qe(1, iNumLoads);
          ColumnVector n1(3), n2(3);
          Matrix dHf_dr(3, 3 * iNumDof), dHf_ds(3, 3 * iNumDof), QA(iNumDof, iNumLoads, 0.);

          for (octave_idx_type i = 0; i < oIntegRule.iGetNumEvalPoints(); ++i) {
               const double alpha = oIntegRule.dGetWeight(i);

               for (octave_idx_type j = 0; j < iNumDir; ++j) {
                    rv.xelem(j) = oIntegRule.dGetPosition(i, j);
               }

               DisplacementInterpMatrixDerR(rv, dHf_dr);
               DisplacementInterpMatrixDerS(rv, dHf_ds);
               ScalarInterpMatrix(rv, HA);

               SurfaceTangentVector(dHf_dr, n1);
               SurfaceTangentVector(dHf_ds, n2);

               const double detJA = JacobianDet(n1, n2);

               info.Add(MeshInfo::JACOBIAN_DET_A, detJA);

               for (octave_idx_type l = 0; l < iNumLoads; ++l) {
                    double HA_qel = 0.;

                    for (octave_idx_type m = 0; m < iNumNodes; ++m) {
                         HA_qel += HA.xelem(m) * qe.xelem(l + iNumLoads * m);
                    }

                    HA_qe.xelem(l) = HA_qel * detJA * alpha;
               }

               for (octave_idx_type l = 0; l < iNumLoads; ++l) {
                    for (octave_idx_type m = 0; m < iNumDof; ++m) {
                         QA.xelem(m + iNumDof * l) += HA.xelem(m) * HA_qe.xelem(l);
                    }
               }
          }

          for (octave_idx_type j = 0; j < iNumLoads; ++j) {
               for (octave_idx_type i = 0; i < iNumDof; ++i) {
                    mat.Insert(QA.xelem(i + iNumDof * j), dofidx.xelem(i), j + colidx);
               }
          }
     }

     virtual octave_idx_type iGetWorkSpaceSize(FemMatrixType eMatType) const override {
          switch (eMatType) {
          case VEC_LOAD_THERMAL:
               return iGetNumDof() * qe.rows();

          default:
               return 0;
          }
     }

     octave_idx_type iGetNumDof() const {
          return nodes.numel();
     }

private:
     const Matrix qe;
     const octave_idx_type colidx;
};

template <typename SHAPE_FUNC, typename BASE>
class SurfaceElementImpl: public BASE {
public:
     template <typename ...Args>
     SurfaceElementImpl(ElementTypes::TypeId eltype, octave_idx_type id, const Matrix& X, const Material* material, const int32NDArray& nodes, const Args&... args)
          :BASE(eltype, id, X, material, nodes, args...) {
          FEM_ASSERT(nodes.numel() == SHAPE_FUNC::iGetNumNodes());
     }

     static void AllocIntegrationRule(Element::FemMatrixType eMatType) {
          SHAPE_FUNC::AllocIntegrationRule(eMatType);
     }

     virtual const IntegrationRule& GetIntegrationRule(Element::FemMatrixType eMatType) const override final {
          return SHAPE_FUNC::GetIntegrationRule(eMatType);
     }

protected:
     virtual void ScalarInterpMatrix(const ColumnVector& rv, Matrix& HA, octave_idx_type irow) const override final {
          SHAPE_FUNC::ScalarInterpMatrix(rv, HA, irow);
     }

     virtual void DisplacementInterpMatrix(const ColumnVector& rv, Matrix& Hf) const override final {
          SHAPE_FUNC::VectorInterpMatrix(rv, Hf);
     }

     virtual void DisplacementInterpMatrixDerR(const ColumnVector& rv, Matrix& dHf_dr) const override final {
          SHAPE_FUNC::VectorInterpMatrixDerR(rv, dHf_dr);
     }

     virtual void DisplacementInterpMatrixDerS(const ColumnVector& rv, Matrix& dHf_ds) const override final {
          SHAPE_FUNC::VectorInterpMatrixDerS(rv, dHf_ds);
     }
};

typedef SurfaceElementImpl<ShapeIso4, PressureLoad> PressureLoadIso4;
typedef SurfaceElementImpl<ShapeQuad8, PressureLoad> PressureLoadQuad8;
typedef SurfaceElementImpl<ShapeQuad8r, PressureLoad> PressureLoadQuad8r;
typedef SurfaceElementImpl<ShapeQuad9, PressureLoad> PressureLoadQuad9;
typedef SurfaceElementImpl<ShapeTria6, PressureLoad> PressureLoadTria6;
typedef SurfaceElementImpl<ShapeTria6H, PressureLoad> PressureLoadTria6H;
typedef SurfaceElementImpl<ShapeTria10, PressureLoad> PressureLoadTria10;

typedef SurfaceElementImpl<ShapeIso4, ThermalConvectionBC> ThermalConvectionBCIso4;
typedef SurfaceElementImpl<ShapeQuad8, ThermalConvectionBC> ThermalConvectionBCQuad8;
typedef SurfaceElementImpl<ShapeQuad8r, ThermalConvectionBC> ThermalConvectionBCQuad8r;
typedef SurfaceElementImpl<ShapeQuad9, ThermalConvectionBC> ThermalConvectionBCQuad9;
typedef SurfaceElementImpl<ShapeTria6, ThermalConvectionBC> ThermalConvectionBCTria6;
typedef SurfaceElementImpl<ShapeTria6H, ThermalConvectionBC> ThermalConvectionBCTria6H;
typedef SurfaceElementImpl<ShapeTria10, ThermalConvectionBC> ThermalConvectionBCTria10;

typedef SurfaceElementImpl<ShapeIso4, HeatSource> HeatSourceIso4;
typedef SurfaceElementImpl<ShapeQuad8, HeatSource> HeatSourceQuad8;
typedef SurfaceElementImpl<ShapeQuad8r, HeatSource> HeatSourceQuad8r;
typedef SurfaceElementImpl<ShapeQuad9, HeatSource> HeatSourceQuad9;
typedef SurfaceElementImpl<ShapeTria6, HeatSource> HeatSourceTria6;
typedef SurfaceElementImpl<ShapeTria6H, HeatSource> HeatSourceTria6H;
typedef SurfaceElementImpl<ShapeTria10, HeatSource> HeatSourceTria10;

typedef SurfaceElementImpl<ShapeIso4, ParticleVelocityBC> ParticleVelocityBCIso4;
typedef SurfaceElementImpl<ShapeQuad8, ParticleVelocityBC> ParticleVelocityBCQuad8;
typedef SurfaceElementImpl<ShapeQuad8r, ParticleVelocityBC> ParticleVelocityBCQuad8r;
typedef SurfaceElementImpl<ShapeQuad9, ParticleVelocityBC> ParticleVelocityBCQuad9;
typedef SurfaceElementImpl<ShapeTria6, ParticleVelocityBC> ParticleVelocityBCTria6;
typedef SurfaceElementImpl<ShapeTria6H, ParticleVelocityBC> ParticleVelocityBCTria6H;
typedef SurfaceElementImpl<ShapeTria10, ParticleVelocityBC> ParticleVelocityBCTria10;

typedef SurfaceElementImpl<ShapeIso4, AcousticImpedanceBC> AcousticImpedanceBCIso4;
typedef SurfaceElementImpl<ShapeQuad8, AcousticImpedanceBC> AcousticImpedanceBCQuad8;
typedef SurfaceElementImpl<ShapeQuad8r, AcousticImpedanceBC> AcousticImpedanceBCQuad8r;
typedef SurfaceElementImpl<ShapeQuad9, AcousticImpedanceBC> AcousticImpedanceBCQuad9;
typedef SurfaceElementImpl<ShapeTria6, AcousticImpedanceBC> AcousticImpedanceBCTria6;
typedef SurfaceElementImpl<ShapeTria6H, AcousticImpedanceBC> AcousticImpedanceBCTria6H;
typedef SurfaceElementImpl<ShapeTria10, AcousticImpedanceBC> AcousticImpedanceBCTria10;

typedef SurfaceElementImpl<ShapeIso4, AcousticBoundary> AcousticBoundaryIso4;
typedef SurfaceElementImpl<ShapeQuad8, AcousticBoundary> AcousticBoundaryQuad8;
typedef SurfaceElementImpl<ShapeQuad8r, AcousticBoundary> AcousticBoundaryQuad8r;
typedef SurfaceElementImpl<ShapeQuad9, AcousticBoundary> AcousticBoundaryQuad9;
typedef SurfaceElementImpl<ShapeTria6, AcousticBoundary> AcousticBoundaryTria6;
typedef SurfaceElementImpl<ShapeTria6H, AcousticBoundary> AcousticBoundaryTria6H;
typedef SurfaceElementImpl<ShapeTria10, AcousticBoundary> AcousticBoundaryTria10;

typedef SurfaceElementImpl<ShapeIso4, FluidStructInteract> FluidStructInteractIso4;
typedef SurfaceElementImpl<ShapeQuad8, FluidStructInteract> FluidStructInteractQuad8;
typedef SurfaceElementImpl<ShapeQuad8r, FluidStructInteract> FluidStructInteractQuad8r;
typedef SurfaceElementImpl<ShapeQuad9, FluidStructInteract> FluidStructInteractQuad9;
typedef SurfaceElementImpl<ShapeTria6, FluidStructInteract> FluidStructInteractTria6;
typedef SurfaceElementImpl<ShapeTria6H, FluidStructInteract> FluidStructInteractTria6H;
typedef SurfaceElementImpl<ShapeTria10, FluidStructInteract> FluidStructInteractTria10;

class StructForce: public Element {
public:
     StructForce(ElementTypes::TypeId eltype, octave_idx_type id, const Matrix& X, const Material* material, const int32NDArray& nodes, octave_idx_type colidx, const Matrix& loads)
          :Element(eltype, id, X, material, nodes),
           loads(loads),
           colidx(colidx) {

          FEM_ASSERT(X.rows() == 3);
          FEM_ASSERT(loads.columns() == 3 || loads.columns() == 6);
          FEM_ASSERT(loads.rows() == nodes.rows());
     }

     StructForce(const StructForce& oElem)
          :Element(oElem.eltype, oElem.id, oElem.X, oElem.material, oElem.nodes), loads(oElem.loads), colidx(oElem.colidx) {
     }

     virtual void Assemble(MatrixAss& mat, MeshInfo& info, const DofMap& dof, FemMatrixType eMatType) const override {
          switch (eMatType) {
          case VEC_LOAD_CONSISTENT:
          case VEC_LOAD_LUMPED:
          case VEC_LOAD_FLUID_STRUCT:
               break;
          default:
               return;
          }

          const octave_idx_type iNumComp = loads.columns();
          const octave_idx_type iNumNodes = loads.rows();

          for (octave_idx_type j = 0; j < iNumComp; ++j) {
               for (octave_idx_type i = 0; i < iNumNodes; ++i) {
                    const octave_idx_type inode = nodes.xelem(i).value() - 1;

                    mat.Insert(loads.xelem(i + iNumNodes * j), dof.GetNodeDofIndex(inode, DofMap::NDOF_DISPLACEMENT, j), colidx);
               }
          }
     }

     virtual octave_idx_type iGetWorkSpaceSize(FemMatrixType eMatType) const override {
          switch (eMatType) {
          case VEC_LOAD_CONSISTENT:
          case VEC_LOAD_LUMPED:
          case VEC_LOAD_FLUID_STRUCT:
               return loads.rows() * loads.columns();
          default:
               return 0;
          }
     }

     static void AllocIntegrationRule(Element::FemMatrixType eMatType) {
     }

private:
     const Matrix loads;
     const octave_idx_type colidx;
};

class ElementBlockBase {
public:
     explicit ElementBlockBase(ElementTypes::TypeId eltype)
          :eltype(eltype) {
     }

     virtual ~ElementBlockBase() {
     }

     virtual octave_idx_type iGetWorkSpaceSize(Element::FemMatrixType eMatType) const=0;
     virtual void Assemble(MatrixAss& oMatAss, MeshInfo& info, const DofMap& oDof, Element::FemMatrixType eMatType, const ParallelOptions& oParaOpt) const=0;
     virtual void PostProcElem(Element::FemMatrixType eMatType, PostProcData& oSolution, const ParallelOptions& oParaOpt) const=0;
     virtual double dGetMass() const=0;
     virtual bool bNeedMatrixInfo(Element::FemMatrixType eMatType) const=0;

     template <typename ElementType, typename... Args>
     inline void Insert(octave_idx_type id, const Matrix& X, const Material* material, const int32NDArray& nodes, const Args&... args);

     ElementTypes::TypeId GetElementType() const { return eltype; }
     virtual void Extract(octave_idx_type& idx, octave_map& sElem) const=0;
     virtual octave_idx_type iGetNumElem() const=0;
     virtual octave_idx_type iGetNumCollocPoints(octave_idx_type id, Element::FemMatrixType eMatType) const=0;

protected:
     const ElementTypes::TypeId eltype;
};

template <typename ElementType>
class ElementBlock: public ElementBlockBase {
public:
     explicit ElementBlock(ElementTypes::TypeId eltype, octave_idx_type iNumElem = 0)
          :ElementBlockBase(eltype) {
          rgElements.reserve(iNumElem);
     }

     template <typename... Args>
     ElementBlock(ElementTypes::TypeId eltype, const int32NDArray& elements, const Matrix& nodes, octave_idx_type inumcoord, const int32NDArray& materials, const vector<Material>& rgMaterials, const Args&... args)
          :ElementBlockBase(eltype) {

          const octave_idx_type iNumElem = elements.rows();
          const octave_idx_type iNumNodesElem = elements.columns();
          const octave_idx_type iNumNodes = nodes.rows();

          rgElements.reserve(iNumElem);

          Matrix X_e(inumcoord, iNumNodesElem);
          int32NDArray nodes_e(dim_vector(iNumNodesElem, 1));

          for (octave_idx_type i = 0; i < iNumElem; ++i) {
               X_e.make_unique();
               nodes_e.make_unique();

               for (octave_idx_type j = 0; j < iNumNodesElem; ++j) {
                    nodes_e.xelem(j) = elements.xelem(i + iNumElem * j);
                    const octave_idx_type jnode = nodes_e.xelem(j).value() - 1;

                    for (octave_idx_type k = 0; k < X_e.rows(); ++k) {
                         X_e.xelem(k + inumcoord * j) = nodes.xelem(jnode + iNumNodes * k);
                    }
               }

               const Material* material = nullptr; // Some elements like RBE3 do not need a material

               if (materials.xelem(i).value() > 0) {
                    material = &rgMaterials[materials.xelem(i).value() - 1];
               }

               rgElements.emplace_back(eltype, i + 1, X_e, material, nodes_e, args...);
          }
     }

     template <typename... Args>
     void Insert(octave_idx_type id, const Matrix& X, const Material* material, const int32NDArray& nodes, const Args&... args) {
          rgElements.emplace_back(eltype, id, X, material, nodes, args...);
     }

     octave_idx_type iGetWorkSpaceSize(Element::FemMatrixType eMatType) const override {
          octave_idx_type iWorkSpace = 0;

          for (auto i = rgElements.begin(); i != rgElements.end(); ++i) {
               const octave_idx_type iCurrWorkSpace = i->ElementType::iGetWorkSpaceSize(eMatType);

               FEM_TRACE("eMatType=" << eMatType << " element=" << i - rgElements.begin() + 1  << " iGetWorkSpaceSize():" << iCurrWorkSpace << "\n");

               FEM_ASSERT(iCurrWorkSpace >= 0);

               iWorkSpace += iCurrWorkSpace;
          }

          return iWorkSpace;
     }

     void Assemble(MatrixAss& oMatAss, MeshInfo& oMeshInfo, const DofMap& oDof, Element::FemMatrixType eMatType, const ParallelOptions& oParaOpt) const override {
          // Allocate integration rule in advance to avoid race condition
          ElementType::AllocIntegrationRule(eMatType);

          const octave_idx_type iNumThreads = oParaOpt.iGetNumThreadsAss();
          const octave_idx_type iWorkSpace = iGetWorkSpaceSize(eMatType);

          if (iNumThreads > 1 && iGetNumElem() >= oParaOpt.iGetMultiThreadThreshold() && iWorkSpace) {
               vector<ThreadData> rgThreadData;

               rgThreadData.reserve(iNumThreads);

               for (octave_idx_type i = 0; i < iNumThreads; ++i) {
                    rgThreadData.emplace_back(oMeshInfo);
               }

               for (const auto& oElem: rgElements) {
                    oElem.ResetElemAssDone();
               }

               const auto pThreadFunc = [this, &oMatAss, &oDof, eMatType] (ThreadData* const pThreadData) {
                                             try {
                                                  for (const auto& oElem: rgElements) {
                                                       OCTAVE_QUIT;

                                                       if (!oElem.bSetElemAssDone()) {
                                                            continue;
                                                       }

                                                       oElem.ElementType::Assemble(oMatAss, pThreadData->oMeshInfo, oDof, eMatType);
                                                  }
                                             } catch (...) {
                                                  pThreadData->pExcept = std::current_exception();
                                             }
                                        };

               vector<std::thread> rgThreads;

               rgThreads.reserve(iNumThreads - 1);

               try {
                    for (octave_idx_type i = 1; i < iNumThreads; ++i) {
                         rgThreads.emplace_back(pThreadFunc, &rgThreadData[i]);
                    }

                    pThreadFunc(&rgThreadData[0]);
               } catch (...) {
                    for (auto& oThread: rgThreads) {
                         oThread.join();
                    }

                    throw;
               }

               for (auto& oThread: rgThreads) {
                    oThread.join();
               }

               for (const ThreadData& oThreadData: rgThreadData) {
                    oMeshInfo.Add(oThreadData.oMeshInfo);

                    if (oThreadData.pExcept) {
                         std::rethrow_exception(oThreadData.pExcept);
                    }
               }
          } else {
               for (const auto& oElem: rgElements) {
                    oElem.ElementType::Assemble(oMatAss, oMeshInfo, oDof, eMatType);

                    OCTAVE_QUIT;
               }
          }
     }

     void PostProcElem(Element::FemMatrixType eMatType, PostProcData& oSolution, const ParallelOptions& oParaOpt) const final override {
          // Allocate integration rule in advance to avoid race condition
          ElementType::AllocIntegrationRule(eMatType);

          const octave_idx_type iNumThreads = oParaOpt.iGetNumThreadsAss();

          if (iNumThreads > 1 && iGetNumElem() >= oParaOpt.iGetMultiThreadThreshold()) {
               vector<std::exception_ptr> rgThreadData;

               rgThreadData.reserve(iNumThreads);

               for (octave_idx_type i = 0; i < iNumThreads; ++i) {
                    rgThreadData.emplace_back(nullptr);
               }

               for (const auto& oElem: rgElements) {
                    oElem.ResetElemAssDone();
               }

               const auto pThreadFunc = [this, &oSolution, eMatType] (std::exception_ptr* const pThreadData) {
                                             try {
                                                  for (const auto& oElem: rgElements) {
                                                       OCTAVE_QUIT;

                                                       if (!oElem.bSetElemAssDone()) {
                                                            continue;
                                                       }

                                                       oElem.ElementType::PostProcElem(eMatType, oSolution);
                                                  }
                                             } catch (...) {
                                                  *pThreadData = std::current_exception();
                                             }
                                        };

               vector<std::thread> rgThreads;

               rgThreads.reserve(iNumThreads - 1);

               try {
                    for (octave_idx_type i = 1; i < iNumThreads; ++i) {
                         rgThreads.emplace_back(pThreadFunc, &rgThreadData[i]);
                    }

                    pThreadFunc(&rgThreadData[0]);
               } catch (...) {
                    for (auto& oThread: rgThreads) {
                         oThread.join();
                    }

                    throw;
               }

               for (auto& oThread: rgThreads) {
                    oThread.join();
               }

               for (const std::exception_ptr& pThreadData: rgThreadData) {
                    if (pThreadData) {
                         std::rethrow_exception(pThreadData);
                    }
               }
          } else {
               for (const auto& oElem: rgElements) {
                    oElem.ElementType::PostProcElem(eMatType, oSolution);

                    OCTAVE_QUIT;
               }
          }
     }

     double dGetMass() const override {
          ElementType::AllocIntegrationRule(Element::MAT_MASS);

          double dm = 0.;

          for (const auto& oElem: rgElements) {
               dm += oElem.ElementType::dGetMass();

               OCTAVE_QUIT;
          }

          return dm;
     }

     virtual bool bNeedMatrixInfo(Element::FemMatrixType eMatType) const override {
          return ElementType::bNeedMatrixInfo(eMatType);
     }

     void Reserve(octave_idx_type iNumElem) {
          rgElements.reserve(iNumElem);
     }

     virtual octave_idx_type iGetNumElem() const override {
          return rgElements.size();
     }

     virtual void Extract(octave_idx_type& idx, octave_map& sElem) const override {
          for (const auto& oElem: rgElements) {
               FEM_ASSERT(sElem.numel() > idx);
               oElem.ElementType::Extract(idx, sElem);
          }
     }

     virtual octave_idx_type iGetNumCollocPoints(octave_idx_type id, Element::FemMatrixType eMatType) const override {
          return rgElements[id].iGetNumCollocPoints(eMatType);
     }
private:
     struct ThreadData {
          explicit ThreadData(const MeshInfo& oMeshInfo)
               :oMeshInfo{oMeshInfo}, pExcept{nullptr} {
          }

          MeshInfo oMeshInfo;
          std::exception_ptr pExcept;
     };

     vector<ElementType> rgElements;
};

template <typename ElementType, typename... Args>
void ElementBlockBase::Insert(octave_idx_type id, const Matrix& X, const Material* material, const int32NDArray& nodes, const Args&... args) {
     typedef ElementBlock<ElementType> ElemBlockType;

     FEM_ASSERT(dynamic_cast<ElemBlockType*>(this) == static_cast<ElemBlockType*>(this));

     static_cast<ElemBlockType*>(this)->Insert(id, X, material, nodes, args...);
}

template <typename PressureElemType>
void InsertPressureElem(const ElementTypes::TypeId eltype, const Matrix& nodes, const octave_map& load_case, const char* const pszElemName, const octave_idx_type iNumNodesElem, vector<std::unique_ptr<ElementBlockBase> >& rgElemBlocks) {
     const octave_idx_type iNumNodes = nodes.rows();
     const auto iter_pressure = load_case.seek("pressure");

     if (iter_pressure != load_case.end()) {
          const Cell cell_pressure = load_case.contents(iter_pressure);

          FEM_ASSERT(cell_pressure.numel() == load_case.numel());

          std::unique_ptr<ElementBlock<PressureElemType> > pElem{nullptr};

          for (octave_idx_type j = 0; j < 2; ++j) {
               octave_idx_type iNumElements = 0;

               for (octave_idx_type i = 0; i < cell_pressure.numel(); ++i) {
                    if (cell_pressure(i).isstruct()) {
                         if (!(cell_pressure(i).numel() == 1)) {
                              throw std::runtime_error("pressure load: pressure must be a scalar struct");
                         }

                         const octave_scalar_map pressure = cell_pressure.xelem(i).scalar_map_value();

                         const auto iter_elem_type = pressure.seek(pszElemName);

                         if (iter_elem_type == pressure.end()) {
                              continue;
                         }

                         const octave_value ov_elem_type = pressure.contents(iter_elem_type);

                         if (!(ov_elem_type.isstruct() && ov_elem_type.numel() == 1)) {
                              throw std::runtime_error("pressure load: invalid entry in load_case.pressure");
                         }

                         const octave_scalar_map elem_type = ov_elem_type.scalar_map_value();

                         const auto iter_elements = elem_type.seek("elements");

                         if (iter_elements == elem_type.end()) {
                              throw std::runtime_error("pressure load: field \"elements\" not found in struct pressure");
                         }

                         const auto iter_p = elem_type.seek("p");

                         if (iter_p == elem_type.end()) {
                              throw std::runtime_error("pressure load: field \"p\" not found in struct pressure");
                         }

                         const octave_value ov_elements = elem_type.contents(iter_elements);

                         if (!(ov_elements.is_matrix_type() && ov_elements.OV_ISINTEGER())) {
                              throw std::runtime_error("pressure load: pressure.elements must be an integer array");
                         }

                         const int32NDArray elements = ov_elements.int32_array_value();

                         if (elements.columns() != iNumNodesElem) {
                              throw std::runtime_error("pressure load: pressure.elements number of columns do not match");
                         }

                         const octave_idx_type iNumElem = elements.rows();

                         const octave_value ov_p = elem_type.contents(iter_p);

                         if (!(ov_p.is_matrix_type() && ov_p.OV_ISREAL())) {
                              throw std::runtime_error("pressure load: pressure.p must be a real matrix");
                         }

                         const Matrix p = ov_p.matrix_value();

                         if (p.columns() != iNumNodesElem || p.rows() != iNumElem) {
                              throw std::runtime_error("pressure load: pressure.p must have the same shape like pressure.elements");
                         }

                         switch (j) {
                         case 0:
                              iNumElements += elements.rows();
                              break;

                         case 1: {
                              Matrix X(3, iNumNodesElem);

                              for (octave_idx_type k = 0; k < iNumElem; ++k) {
                                   X.make_unique();

                                   for (octave_idx_type l = 0; l < iNumNodesElem; ++l) {
                                        for (octave_idx_type m = 0; m < 3; ++m) {
                                             const octave_idx_type inode = elements.xelem(k + iNumElem * l).value() - 1;

                                             if (inode < 0 || inode >= nodes.rows()) {
                                                  throw std::runtime_error("pressure load: node index out of range in pressure.elements");
                                             }

                                             X.xelem(m + 3 * l) = nodes.xelem(inode + iNumNodes * m);
                                        }
                                   }

                                   pElem->Insert(++iNumElements,
                                                 X,
                                                 nullptr,
                                                 elements.index(idx_vector::make_range(k, 1, 1),
                                                                idx_vector::make_range(0, 1, elements.columns())),
                                                 p.row(k),
                                                 i + 1);
                              }
                         } break;

                         default:
                              FEM_ASSERT(false);
                         }
                    }
               }

               if (iNumElements) {
                    if (j == 0) {
                         pElem.reset(new ElementBlock<PressureElemType>(eltype));
                         pElem->Reserve(iNumElements);
                    }
               } else {
                    break;
               }
          }

          if (pElem) {
               rgElemBlocks.emplace_back(std::move(pElem));
          }
     }
}

template <typename ConvectionElemType>
void InsertThermalConvElem(const ElementTypes::TypeId eltype, const Matrix& nodes, const octave_scalar_map& elements, const octave_map& load_case, const char* const pszElemName, const octave_idx_type iNumNodesElem, vector<std::unique_ptr<ElementBlockBase> >& rgElemBlocks) {
     const octave_idx_type iNumNodes = nodes.rows();
     const octave_idx_type iNumLoads = load_case.numel();

     const auto iter_convection = elements.seek("convection");

     if (iter_convection == elements.end()) {
          return;
     }

     const octave_value ov_convection = elements.contents(iter_convection);

     if (!(ov_convection.isstruct() && ov_convection.numel() == 1)) {
          throw std::runtime_error("thermal convection: mesh.elements.convection must be a scalar struct");
     }

     const octave_scalar_map m_convection = ov_convection.scalar_map_value();

     const auto iter_elem_type = m_convection.seek(pszElemName);

     if (iter_elem_type == m_convection.end()) {
          return;
     }

     const octave_value ov_elem_type = m_convection.contents(iter_elem_type);

     if (!(ov_elem_type.isstruct() && ov_elem_type.numel() == 1)) {
          throw std::runtime_error("thermal convection: mesh.elements.convection."s + pszElemName + " must be a scalar struct");
     }

     const octave_scalar_map m_elem_type = ov_elem_type.scalar_map_value();

     const auto iter_elnodes = m_elem_type.seek("nodes");

     if (iter_elnodes == m_elem_type.end()) {
          throw std::runtime_error("thermal convection: missing field mesh.elements.convection."s + pszElemName + ".nodes");
     }

     const octave_value ov_elnodes = m_elem_type.contents(iter_elnodes);

     if (!(ov_elnodes.is_matrix_type() && ov_elnodes.isinteger() && ov_elnodes.columns() == iNumNodesElem)) {
          throw std::runtime_error("thermal convection: mesh.elements.convection."s + pszElemName + ".nodes must be an integer matrix");
     }

     const int32NDArray elnodes = ov_elnodes.int32_array_value();
     const octave_idx_type iNumElem = elnodes.rows();

     const auto iter_h = m_elem_type.seek("h");

     if (iter_h == m_elem_type.end()) {
          throw std::runtime_error("thermal convection: missing field mesh.elements.convection."s + pszElemName + ".h");
     }

     const octave_value ov_h = m_elem_type.contents(iter_h);

     if (!(ov_h.is_matrix_type() && ov_h.isreal() && ov_h.rows() == iNumElem && ov_h.columns() == iNumNodesElem)) {
          throw std::runtime_error("thermal convection: mesh.elements.convection."s + pszElemName
                                   + ".h must be a real matrix with the same size like mesh.elements.convection."
                                   + pszElemName + ".nodes");
     }

     const Matrix h = ov_h.matrix_value();

     NDArray theta(dim_vector(iNumLoads, iNumNodesElem, iNumElem), 0.);

     const auto iter_conv_load = load_case.seek("convection");

     if (iter_conv_load != load_case.end()) {
          const Cell cell_conv_load = load_case.contents(iter_conv_load);

          FEM_ASSERT(cell_conv_load.numel() == theta.rows());

          for (octave_idx_type k = 0; k < iNumLoads; ++k) {
               const octave_value ov_conv_load = cell_conv_load.xelem(k);

               if (!(ov_conv_load.isstruct() && ov_conv_load.numel() == 1)) {
                    throw std::runtime_error("thermal convection: load_case.convection must be a scalar struct");
               }

               const octave_scalar_map m_conv_load = ov_conv_load.scalar_map_value();

               const auto iter_conv_load_elem = m_conv_load.seek(pszElemName);

               if (iter_conv_load_elem == m_conv_load.end()) {
                    continue;
               }

               const octave_value ov_conv_load_elem = m_conv_load.contents(iter_conv_load_elem);

               if (!(ov_conv_load_elem.isstruct() && ov_conv_load_elem.numel() == 1)) {
                    throw std::runtime_error("thermal convection: load_case.convection."s + pszElemName + " must be a scalar struct");
               }

               const octave_scalar_map m_conv_load_elem = ov_conv_load_elem.scalar_map_value();

               const auto iter_theta = m_conv_load_elem.seek("theta");

               if (iter_theta == m_conv_load_elem.end()) {
                    throw std::runtime_error("thermal convection: missing field load_case.convection."s + pszElemName + ".theta");
               }

               const octave_value ov_thetak = m_conv_load_elem.contents(iter_theta);

               if (!(ov_thetak.is_matrix_type() && ov_thetak.isreal() && ov_thetak.rows() == iNumElem && ov_thetak.columns() == iNumNodesElem)) {
                    throw std::runtime_error("thermal convection: load_case.convection."s + pszElemName + ".theta must be a real matrix with the same dimensions like mesh.elements.convection."s + pszElemName + ".nodes");
               }

               const Matrix thetak = ov_thetak.matrix_value();

               for (octave_idx_type j = 0; j < iNumElem; ++j) {
                    for (octave_idx_type i = 0; i < iNumNodesElem; ++i) {
                         theta.xelem(k + iNumLoads * (i + iNumNodesElem * j)) = thetak.xelem(j + iNumElem * i);
                    }
               }
          }
     }

     NDArray X(dim_vector(3, iNumNodesElem, iNumElem));

     for (octave_idx_type k = 0; k < iNumElem; ++k) {
          for (octave_idx_type l = 0; l < iNumNodesElem; ++l) {
               const octave_idx_type inode = elnodes.xelem(k + iNumElem * l).value() - 1;

               if (inode < 0 || inode >= nodes.rows()) {
                    throw std::runtime_error("thermal convection: node index out of range in mesh.elements.convection."s + pszElemName + ".nodes");
               }

               for (octave_idx_type m = 0; m < 3; ++m) {
                    X.xelem(m + 3 * (l + iNumNodesElem * k)) = nodes.xelem(inode + iNumNodes * m);
               }
          }
     }

     std::unique_ptr<ElementBlock<ConvectionElemType> > pElem{new ElementBlock<ConvectionElemType>(eltype)};

     pElem->Reserve(elnodes.rows());

     for (octave_idx_type k = 0; k < iNumElem; ++k) {
          const Matrix Xk = X.linear_slice(3 * iNumNodesElem * k, 3 * iNumNodesElem * (k + 1)).reshape(dim_vector(3, iNumNodesElem));
          const Matrix thetak = theta.linear_slice(iNumLoads * iNumNodesElem * k, iNumLoads * iNumNodesElem * (k + 1)).reshape(dim_vector(iNumLoads, iNumNodesElem));

          pElem->Insert(k + 1,
                        Xk,
                        nullptr,
                        elnodes.index(idx_vector::make_range(k, 1, 1),
                                      idx_vector::make_range(0, 1, iNumNodesElem)),
                        thetak,
                        h.row(k));
     }

     rgElemBlocks.emplace_back(std::move(pElem));
}

template <typename SourceElemType>
void InsertHeatSourceElem(const ElementTypes::TypeId eltype, const Matrix& nodes, const octave_map& load_case, const char* const pszElemName, const octave_idx_type iNumNodesElem, vector<std::unique_ptr<ElementBlockBase> >& rgElemBlocks) {
     const octave_idx_type iNumNodes = nodes.rows();
     const auto iter_heat_source = load_case.seek("heat_source");

     if (iter_heat_source == load_case.end()) {
          return;
     }

     const Cell cell_heat_source = load_case.contents(iter_heat_source);

     std::unique_ptr<ElementBlock<SourceElemType> > pElem;

     octave_idx_type iNumElemTot = 0;

     for (octave_idx_type i = 0; i < 2; ++i) {
          for (octave_idx_type j = 0; j < cell_heat_source.numel(); ++j) {
               const octave_value ov_heat_source = cell_heat_source.xelem(j);

               if (!(ov_heat_source.isstruct() && ov_heat_source.numel() == 1)) {
                    throw std::runtime_error("heat source: load_case.heat_source must be a scalar struct");
               }

               const octave_scalar_map m_heat_source = ov_heat_source.scalar_map_value();

               const auto iter_heat_source_elem = m_heat_source.seek(pszElemName);

               if (iter_heat_source_elem == m_heat_source.end()) {
                    continue;
               }

               const octave_value ov_heat_source_elem = m_heat_source.contents(iter_heat_source_elem);

               if (!(ov_heat_source_elem.isstruct() && ov_heat_source_elem.numel() == 1)) {
                    throw std::runtime_error("heat source: load_case.heat_source."s + pszElemName + " must be a scalar struct");
               }

               const octave_scalar_map m_heat_source_elem = ov_heat_source_elem.scalar_map_value();

               const auto iter_elnodes = m_heat_source_elem.seek("nodes");

               if (iter_elnodes == m_heat_source_elem.end()) {
                    throw std::runtime_error("heat source: missing field load_case.heat_source."s + pszElemName + ".nodes");
               }

               const octave_value ov_elnodes = m_heat_source_elem.contents(iter_elnodes);

               if (!(ov_elnodes.is_matrix_type() && ov_elnodes.isinteger() && ov_elnodes.columns() == iNumNodesElem)) {
                    throw std::runtime_error("heat source: load_case.heat_source."s + pszElemName + ".nodes must be an integer matrix");
               }

               const int32NDArray elnodes = ov_elnodes.int32_array_value();

               const octave_idx_type iNumElem = elnodes.rows();

               if (i == 0) {
                    iNumElemTot += iNumElem;
                    continue;
               }

               const auto iter_q = m_heat_source_elem.seek("q");

               if (iter_q == m_heat_source_elem.end()) {
                    throw std::runtime_error("heat source: missing field load_case.heat_source."s + pszElemName + ".q");
               }

               const octave_value ov_q = m_heat_source_elem.contents(iter_q);

               if (!(ov_q.is_matrix_type() && ov_q.isreal() && ov_q.rows() == iNumElem && ov_q.columns() == iNumNodesElem)) {
                    throw std::runtime_error("heat source: load_case.heat_source."s + pszElemName + ".q must be a real matrix with the same dimensions like load_case.heat_source."s + pszElemName + ".nodes");
               }

               const Matrix q = ov_q.matrix_value();

               NDArray X(dim_vector(3, iNumNodesElem, iNumElem));

               for (octave_idx_type k = 0; k < iNumElem; ++k) {
                    for (octave_idx_type l = 0; l < iNumNodesElem; ++l) {
                         const octave_idx_type inode = elnodes.xelem(k + iNumElem * l).value() - 1;

                         if (inode < 0 || inode >= iNumNodes) {
                              throw std::runtime_error("heat source: node index out of range in load_case.heat_source."s + pszElemName + ".nodes");
                         }

                         for (octave_idx_type m = 0; m < 3; ++m) {
                              X.xelem(m + 3 * (l + iNumNodesElem * k)) = nodes.xelem(inode + iNumNodes * m);
                         }
                    }
               }

               for (octave_idx_type k = 0; k < iNumElem; ++k) {
                    const Matrix Xk = X.linear_slice(3 * iNumNodesElem * k, 3 * iNumNodesElem * (k + 1)).reshape(dim_vector(3, iNumNodesElem));

                    pElem->Insert(k + 1,
                                  Xk,
                                  nullptr,
                                  elnodes.index(idx_vector::make_range(k, 1, 1),
                                                idx_vector::make_range(0, 1, iNumNodesElem)),
                                  q.row(k),
                                  j + 1);
               }
          }

          if (i == 0) {
               pElem.reset(new ElementBlock<SourceElemType>(eltype));
               pElem->Reserve(iNumElemTot);
          }
     }

     rgElemBlocks.emplace_back(std::move(pElem));
}

template <typename VelocityElemType>
void InsertParticleVelocityBC(const ElementTypes::TypeId eltype, const Matrix& nodes, const octave_scalar_map& elements, const std::vector<Material>& rgMaterials, const octave_scalar_map& materials, const octave_map& load_case, const char* const pszElemName, const octave_idx_type iNumNodesElem, vector<std::unique_ptr<ElementBlockBase> >& rgElemBlocks, const DofMap::DomainType eDomain) {
     const octave_idx_type iNumNodes = nodes.rows();
     const octave_idx_type iNumLoads = load_case.numel();
     const auto iter_velocity = elements.seek("particle_velocity");

     if (iter_velocity == elements.end()) {
          return;
     }

     const octave_value ov_velocity = elements.contents(iter_velocity);

     if (!(ov_velocity.isstruct() && ov_velocity.numel() == 1)) {
          throw std::runtime_error("particle velocity: mesh.elements.particle_velocity must be a scalar struct");
     }

     const octave_scalar_map m_velocity = ov_velocity.scalar_map_value();

     const auto iter_elem_type = m_velocity.seek(pszElemName);

     if (iter_elem_type == m_velocity.end()) {
          return;
     }

     const octave_value ov_elem_type = m_velocity.contents(iter_elem_type);

     if (!(ov_elem_type.isstruct() && ov_elem_type.numel() == 1)) {
          throw std::runtime_error("particle velocity: mesh.elements.particle_velocity."s + pszElemName + " must be a scalar struct");
     }

     const octave_scalar_map m_elem_type = ov_elem_type.scalar_map_value();

     const auto iter_elnodes = m_elem_type.seek("nodes");

     if (iter_elnodes == m_elem_type.end()) {
          throw std::runtime_error("particle velocity: missing field mesh.elements.particle_velocity."s + pszElemName + ".nodes");
     }

     const octave_value ov_elnodes = m_elem_type.contents(iter_elnodes);

     if (!(ov_elnodes.is_matrix_type() && ov_elnodes.isinteger() && ov_elnodes.columns() == iNumNodesElem)) {
          throw std::runtime_error("particle velocity: mesh.elements.particle_velocity."s + pszElemName + ".nodes must be an integer matrix");
     }

     const int32NDArray elnodes = ov_elnodes.int32_array_value();
     const octave_idx_type iNumElem = elnodes.rows();

     NDArray X(dim_vector(3, iNumNodesElem, iNumElem));

     for (octave_idx_type k = 0; k < iNumElem; ++k) {
          for (octave_idx_type l = 0; l < iNumNodesElem; ++l) {
               const octave_idx_type inode = elnodes.xelem(k + iNumElem * l).value() - 1;

               if (inode < 0 || inode >= iNumNodes) {
                    throw std::runtime_error("particle velocity: node index out of range in mesh.elements.particle_velocity."s + pszElemName + ".nodes");
               }

               for (octave_idx_type m = 0; m < 3; ++m) {
                    X.xelem(m + 3 * (l + iNumNodesElem * k)) = nodes.xelem(inode + iNumNodes * m);
               }
          }
     }

     const auto iter_vel_mat = materials.seek("particle_velocity");

     if (iter_vel_mat == materials.end()) {
          throw std::runtime_error("mesh.materials.particle_velocity is not defined");
     }

     const octave_scalar_map m_vel_mat = materials.contents(iter_vel_mat).scalar_map_value();

     const auto iter_elem_type_mat = m_vel_mat.seek(pszElemName);

     if (iter_elem_type_mat == m_vel_mat.end()) {
          throw std::runtime_error("particle velocity: mesh.materials.particle_velocity."s + pszElemName + " is not defined");
     }

     const int32NDArray elem_mat = m_vel_mat.contents(iter_elem_type_mat).int32_array_value();

     if (elem_mat.numel() != iNumElem) {
          throw std::runtime_error("particle velocity: invalid number of rows for matrix mesh.materials.particle_velocity."s + pszElemName + " in argument mesh");
     }

     const octave_idx_type inum_materials = rgMaterials.size();

     for (octave_idx_type i = 0; i < iNumElem; ++i) {
          const octave_idx_type imaterial = elem_mat.xelem(i);

          if (imaterial <= 0 || imaterial > inum_materials) {
               throw std::runtime_error("particle velocity: invalid index in matrix mesh.materials.particle_velocity."s + pszElemName + " in argument mesh");
          }
     }

     NDArray vel(dim_vector(iNumLoads, iNumNodesElem, iNumElem), 0.);

     const auto iter_vel_load = load_case.seek("particle_velocity");

     if (iter_vel_load != load_case.end()) {
          const Cell cell_vel_load = load_case.contents(iter_vel_load);

          FEM_ASSERT(cell_vel_load.numel() == iNumLoads);

          for (octave_idx_type k = 0; k < iNumLoads; ++k) {
               const octave_value ov_vel_load = cell_vel_load.xelem(k);

               if (!(ov_vel_load.isstruct() && ov_vel_load.numel() == 1)) {
                    throw std::runtime_error("particle velocity: load_case.particle_velocity must be a scalar struct");
               }

               const octave_scalar_map m_vel_load = ov_vel_load.scalar_map_value();

               const auto iter_vel_load_elem = m_vel_load.seek(pszElemName);

               if (iter_vel_load_elem == m_vel_load.end()) {
                    continue;
               }

               const octave_value ov_vel_load_elem = m_vel_load.contents(iter_vel_load_elem);

               if (!(ov_vel_load_elem.isstruct() && ov_vel_load_elem.numel() == 1)) {
                    throw std::runtime_error("particle velocity: load_case.particle_velocity."s + pszElemName + " must be a scalar struct");
               }

               const octave_scalar_map m_vel_load_elem = ov_vel_load_elem.scalar_map_value();

               const auto iter_vel = m_vel_load_elem.seek("vn");

               if (iter_vel == m_vel_load_elem.end()) {
                    throw std::runtime_error("particle velocity: missing field load_case.particle_velocity."s + pszElemName + ".vn");
               }

               const octave_value ov_vnk = m_vel_load_elem.contents(iter_vel);

               if (!(ov_vnk.is_matrix_type() && ov_vnk.isreal() && ov_vnk.rows() == iNumElem && ov_vnk.columns() == iNumNodesElem)) {
                    throw std::runtime_error("load_case.particle_velocity."s + pszElemName + ".vn must be a real matrix with the same dimensions like mesh.elements.particle_velocity."s + pszElemName + ".nodes");
               }

               const Matrix vk = ov_vnk.matrix_value();

               for (octave_idx_type j = 0; j < iNumElem; ++j) {
                    for (octave_idx_type i = 0; i < iNumNodesElem; ++i) {
                         vel.xelem(k + iNumLoads * (i + iNumNodesElem * j)) = vk.xelem(j + iNumElem * i);
                    }
               }
          }
     }

     std::unique_ptr<ElementBlock<VelocityElemType> > pElem{new ElementBlock<VelocityElemType>(eltype)};

     pElem->Reserve(iNumElem);

     const RowVector coef(iNumNodesElem, eDomain == DofMap::DO_ACOUSTICS ? 1. : -1.);

     for (octave_idx_type k = 0; k < iNumElem; ++k) {
          const Matrix Xk = X.linear_slice(3 * iNumNodesElem * k, 3 * iNumNodesElem * (k + 1)).reshape(dim_vector(3, iNumNodesElem));
          const Matrix velk = vel.linear_slice(iNumLoads * iNumNodesElem * k, iNumLoads * iNumNodesElem * (k + 1)).reshape(dim_vector(iNumLoads, iNumNodesElem));

          FEM_ASSERT(static_cast<size_t>(elem_mat.xelem(k).value() - 1) < rgMaterials.size());

          const Material* const materialk = &rgMaterials[elem_mat.xelem(k).value() - 1];

          pElem->Insert(k + 1,
                        Xk,
                        materialk,
                        elnodes.index(idx_vector::make_range(k, 1, 1),
                                      idx_vector::make_range(0, 1, iNumNodesElem)),
                        velk,
                        coef);
     }

     rgElemBlocks.emplace_back(std::move(pElem));
}

template <typename ImpedanceElemType>
void InsertAcousticImpedanceBC(const ElementTypes::TypeId eltype, const Matrix& nodes, const octave_scalar_map& elements, const std::vector<Material>& rgMaterials, const octave_scalar_map& materials, const char* const pszElemName, const octave_idx_type iNumNodesElem, vector<std::unique_ptr<ElementBlockBase> >& rgElemBlocks, const DofMap::DomainType eDomain) {
     const octave_idx_type iNumNodes = nodes.rows();
     const auto iter_impedance = elements.seek("acoustic_impedance");

     if (iter_impedance == elements.end()) {
          return;
     }

     const octave_value ov_impedance = elements.contents(iter_impedance);

     if (!(ov_impedance.isstruct() && ov_impedance.numel() == 1)) {
          throw std::runtime_error("acoustic impedance: mesh.elements.acoustic_impedance must be a scalar struct");
     }

     const octave_scalar_map m_impedance = ov_impedance.scalar_map_value();

     const auto iter_elem_type = m_impedance.seek(pszElemName);

     if (iter_elem_type == m_impedance.end()) {
          return;
     }

     const octave_value ov_elem_type = m_impedance.contents(iter_elem_type);

     if (!(ov_elem_type.isstruct() && ov_elem_type.numel() == 1)) {
          throw std::runtime_error("acoustic impedance: mesh.elements.acoustic_impedance."s + pszElemName + " must be a scalar struct");
     }

     const octave_scalar_map m_elem_type = ov_elem_type.scalar_map_value();

     const auto iter_elnodes = m_elem_type.seek("nodes");

     if (iter_elnodes == m_elem_type.end()) {
          throw std::runtime_error("acoustic impedance: missing field mesh.elements.acoustic_impedance."s + pszElemName + ".nodes");
     }

     const octave_value ov_elnodes = m_elem_type.contents(iter_elnodes);

     if (!(ov_elnodes.is_matrix_type() && ov_elnodes.isinteger() && ov_elnodes.columns() == iNumNodesElem)) {
          throw std::runtime_error("acoustic impedance: mesh.elements.acoustic_impedance."s + pszElemName + ".nodes must be an integer matrix");
     }

     const int32NDArray elnodes = ov_elnodes.int32_array_value();
     const octave_idx_type iNumElem = elnodes.rows();

     NDArray X(dim_vector(3, iNumNodesElem, iNumElem));

     for (octave_idx_type k = 0; k < iNumElem; ++k) {
          for (octave_idx_type l = 0; l < iNumNodesElem; ++l) {
               octave_idx_type inode = elnodes.xelem(k + iNumElem * l).value() - 1;

               if (inode < 0 || inode >= nodes.rows()) {
                    throw std::runtime_error("acoustic impedance: node index out of range in mesh.elements.acoustic_impedance."s + pszElemName + ".nodes");
               }

               for (octave_idx_type m = 0; m < 3; ++m) {
                    X.xelem(m + 3 * (l + iNumNodesElem * k)) = nodes.xelem(inode + iNumNodes * m);
               }
          }
     }

     const auto iter_z = m_elem_type.seek("z");

     if (iter_z == m_elem_type.end()) {
          throw std::runtime_error("acoustic impedance: missing field mesh.elements.acoustic_impedance."s + pszElemName + ".z");
     }

     const octave_value ov_z = m_elem_type.contents(iter_z);

     if (!(ov_z.is_matrix_type() && ov_z.rows() == iNumElem && ov_z.columns() == iNumNodesElem)) {
          throw std::runtime_error("acoustic impedance: mesh.elements.acoustic_impedance."s + pszElemName + ".z must be a real matrix of the same size "
                                   "like mesh.elements.acoustic_impedance." + pszElemName + ".nodes");
     }

     const ComplexMatrix z = ov_z.complex_matrix_value();

     const auto iter_impe_mat = materials.seek("acoustic_impedance");

     if (iter_impe_mat == materials.end()) {
          throw std::runtime_error("acoustic impedance: mesh.materials.acoustic_impedance is not defined");
     }

     const octave_scalar_map m_impe_mat = materials.contents(iter_impe_mat).scalar_map_value();

     const auto iter_elem_type_mat = m_impe_mat.seek(pszElemName);

     if (iter_elem_type_mat == m_impe_mat.end()) {
          throw std::runtime_error("acoustic impedance: mesh.materials.acoustic_impedance."s + pszElemName + " is not defined");
     }

     const int32NDArray elem_mat = m_impe_mat.contents(iter_elem_type_mat).int32_array_value();

     if (elem_mat.numel() != iNumElem) {
          throw std::runtime_error("acoustic impedance: invalid number of rows for matrix mesh.materials.acoustic_impedance."s + pszElemName + " in argument mesh");
     }

     const octave_idx_type inum_materials = rgMaterials.size();

     for (octave_idx_type i = 0; i < iNumElem; ++i) {
          const octave_idx_type imaterial = elem_mat.xelem(i);

          if (imaterial <= 0 || imaterial > inum_materials) {
               throw std::runtime_error("acoustic impedance: invalid index in matrix mesh.materials.acoustic_impedance."s + pszElemName + " in argument mesh");
          }
     }

     std::unique_ptr<ElementBlock<ImpedanceElemType> > pElem{new ElementBlock<ImpedanceElemType>(eltype)};

     pElem->Reserve(iNumElem);

     const double coef = eDomain == DofMap::DO_ACOUSTICS ? 1. : -1.;

     for (octave_idx_type k = 0; k < iNumElem; ++k) {
          const Matrix Xk = X.linear_slice(3 * iNumNodesElem * k, 3 * iNumNodesElem * (k + 1)).reshape(dim_vector(3, iNumNodesElem));

          FEM_ASSERT(static_cast<size_t>(elem_mat.xelem(k).value() - 1) < rgMaterials.size());

          const Material* const materialk = &rgMaterials[elem_mat.xelem(k).value() - 1];
          RowVector rec_z_re(iNumNodesElem), rec_z_im(iNumNodesElem);

          for (octave_idx_type i = 0; i < iNumNodesElem; ++i) {
               const std::complex<double> rec_zki = coef / z.xelem(k + iNumElem * i);
               rec_z_re.xelem(i) = std::real(rec_zki);
               rec_z_im.xelem(i) = std::imag(rec_zki);
          }

          pElem->Insert(k + 1,
                        Xk,
                        materialk,
                        elnodes.index(idx_vector::make_range(k, 1, 1),
                                      idx_vector::make_range(0, 1, iNumNodesElem)),
                        rec_z_re,
                        rec_z_im);
     }

     rgElemBlocks.emplace_back(std::move(pElem));
}

template <typename BoundaryElemType>
void InsertAcousticBoundary(const ElementTypes::TypeId eltype, const Matrix& nodes, const octave_scalar_map& elements, const std::vector<Material>& rgMaterials, const octave_scalar_map& materials, const char* const pszElemName, const octave_idx_type iNumNodesElem, vector<std::unique_ptr<ElementBlockBase> >& rgElemBlocks) {
     const octave_idx_type iNumNodes = nodes.rows();
     const auto iter_boundary = elements.seek("acoustic_boundary");

     if (iter_boundary == elements.end()) {
          return;
     }

     const octave_value ov_boundary = elements.contents(iter_boundary);

     if (!(ov_boundary.isstruct() && ov_boundary.numel() == 1)) {
          throw std::runtime_error("acoustic boundary: mesh.elements.acoustic_boundary must be a scalar struct");
     }

     const octave_scalar_map m_boundary = ov_boundary.scalar_map_value();

     const auto iter_elem_type = m_boundary.seek(pszElemName);

     if (iter_elem_type == m_boundary.end()) {
          return;
     }

     const octave_value ov_elnodes = m_boundary.contents(iter_elem_type);

     if (!(ov_elnodes.is_matrix_type() && ov_elnodes.isinteger() && ov_elnodes.columns() == iNumNodesElem)) {
          throw std::runtime_error("acoustic boundary: mesh.elements.acoustic_boundary."s + pszElemName + " must be an integer matrix");
     }

     const int32NDArray elnodes = ov_elnodes.int32_array_value();
     const octave_idx_type iNumElem = elnodes.rows();

     NDArray X(dim_vector(3, iNumNodesElem, iNumElem));

     for (octave_idx_type k = 0; k < iNumElem; ++k) {
          for (octave_idx_type l = 0; l < iNumNodesElem; ++l) {
               const octave_idx_type inode = elnodes.xelem(k + iNumElem * l).value() - 1;

               if (inode < 0 || inode >= nodes.rows()) {
                    throw std::runtime_error("acoustic boundary: node index out of range in mesh.elements.acoustic_boundary."s + pszElemName);
               }

               for (octave_idx_type m = 0; m < 3; ++m) {
                    X.xelem(m + 3 * (l + iNumNodesElem * k)) = nodes.xelem(inode + iNumNodes * m);
               }
          }
     }

     const auto iter_bnd_mat = materials.seek("acoustic_boundary");

     if (iter_bnd_mat == materials.end()) {
          throw std::runtime_error("acoustic boundary: mesh.materials.acoustic_boundary is not defined");
     }

     const octave_scalar_map m_bnd_mat = materials.contents(iter_bnd_mat).scalar_map_value();

     const auto iter_elem_type_mat = m_bnd_mat.seek(pszElemName);

     if (iter_elem_type_mat == m_bnd_mat.end()) {
          throw std::runtime_error("acoustic boundary: mesh.materials.acoustic_boundary."s + pszElemName + " is not defined");
     }

     const int32NDArray elem_mat = m_bnd_mat.contents(iter_elem_type_mat).int32_array_value();

     if (elem_mat.numel() != iNumElem) {
          throw std::runtime_error("acoustic boundary: invalid number of rows for matrix mesh.materials.acoustic_boundary."s + pszElemName + " in argument mesh");
     }

     const octave_idx_type inum_materials = rgMaterials.size();

     for (octave_idx_type i = 0; i < elem_mat.numel(); ++i) {
          const octave_idx_type imaterial = elem_mat.xelem(i);

          if (imaterial <= 0 || imaterial > inum_materials) {
               throw std::runtime_error("acoustic boundary: invalid index in matrix mesh.materials.acoustic_boundary."s + pszElemName + " in argument mesh");
          }
     }

     std::unique_ptr<ElementBlock<BoundaryElemType> > pElem{new ElementBlock<BoundaryElemType>(eltype)};

     pElem->Reserve(iNumElem);

     for (octave_idx_type k = 0; k < iNumElem; ++k) {
          const Matrix Xk = X.linear_slice(3 * iNumNodesElem * k, 3 * iNumNodesElem * (k + 1)).reshape(dim_vector(3, iNumNodesElem));

          FEM_ASSERT(static_cast<size_t>(elem_mat.xelem(k).value() - 1) < rgMaterials.size());

          const Material* const materialk = &rgMaterials[elem_mat.xelem(k).value() - 1];

          pElem->Insert(k + 1,
                        Xk,
                        materialk,
                        elnodes.index(idx_vector::make_range(k, 1, 1),
                                      idx_vector::make_range(0, 1, iNumNodesElem)));
     }

     rgElemBlocks.emplace_back(std::move(pElem));
}

template <typename FluidStructElemType>
void InsertFluidStructElem(const ElementTypes::TypeId eltype, const Matrix& nodes, const octave_scalar_map& elements, const char* const pszElemName, const octave_idx_type iNumNodesElem, vector<std::unique_ptr<ElementBlockBase> >& rgElemBlocks) {
     const octave_idx_type iNumNodes = nodes.rows();
     const auto iter_fluid_struct = elements.seek("fluid_struct_interface");

     if (iter_fluid_struct == elements.end()) {
          return;
     }

     const octave_value ov_fluid_struct = elements.contents(iter_fluid_struct);

     if (!(ov_fluid_struct.isstruct() && ov_fluid_struct.numel() == 1)) {
          throw std::runtime_error("fluid struct interface: mesh.elements.fluid_struct_interface must be a scalar struct");
     }

     const octave_scalar_map m_fluid_struct = ov_fluid_struct.scalar_map_value();

     const auto iter_elem_type = m_fluid_struct.seek(pszElemName);

     if (iter_elem_type == m_fluid_struct.end()) {
          return;
     }

     const octave_value ov_elnodes = m_fluid_struct.contents(iter_elem_type);

     if (!(ov_elnodes.is_matrix_type() && ov_elnodes.isinteger() && ov_elnodes.columns() == iNumNodesElem)) {
          throw std::runtime_error("fluid struct interface: mesh.elements.fluid_struct_interface."s + pszElemName + " must be an integer matrix");
     }

     const int32NDArray elnodes = ov_elnodes.int32_array_value();
     const octave_idx_type iNumElem = elnodes.rows();

     NDArray X(dim_vector(3, iNumNodesElem, iNumElem));

     for (octave_idx_type k = 0; k < iNumElem; ++k) {
          for (octave_idx_type l = 0; l < iNumNodesElem; ++l) {
               const octave_idx_type inode = elnodes.xelem(k + iNumElem * l).value() - 1;

               if (inode < 0 || inode >= nodes.rows()) {
                    throw std::runtime_error("fluid struct interface: node index out of range in mesh.elements.fluid_struct_interface."s + pszElemName + ".nodes");
               }

               for (octave_idx_type m = 0; m < 3; ++m) {
                    X.xelem(m + 3 * (l + iNumNodesElem * k)) = nodes.xelem(inode + iNumNodes * m);
               }
          }
     }

     std::unique_ptr<ElementBlock<FluidStructElemType> > pElem{new ElementBlock<FluidStructElemType>{eltype}};

     pElem->Reserve(iNumElem);

     for (octave_idx_type k = 0; k < iNumElem; ++k) {
          const Matrix Xk = X.linear_slice(3 * iNumNodesElem * k, 3 * iNumNodesElem * (k + 1)).reshape(dim_vector(3, iNumNodesElem));

          pElem->Insert(k + 1,
                        Xk,
                        nullptr,
                        elnodes.index(idx_vector::make_range(k, 1, 1),
                                      idx_vector::make_range(0, 1, iNumNodesElem)));
     }

     rgElemBlocks.emplace_back(std::move(pElem));
}

namespace shape_func_util {
     template <typename T>
     struct SelectElemPerSlaveNode {
     };

     template <>
     struct SelectElemPerSlaveNode<ShapeIso4> {
          static constexpr octave_idx_type iNumElemPerSlaveNode = 9;
     };

     template <>
     struct SelectElemPerSlaveNode<ShapeQuad8> {
          static constexpr octave_idx_type iNumElemPerSlaveNode = 9;
     };

     template <>
     struct SelectElemPerSlaveNode<ShapeQuad8r> {
          static constexpr octave_idx_type iNumElemPerSlaveNode = 9;
     };

     template <>
     struct SelectElemPerSlaveNode<ShapeQuad9> {
          static constexpr octave_idx_type iNumElemPerSlaveNode = 9;
     };

     template <>
     struct SelectElemPerSlaveNode<ShapeTria6> {
          static constexpr octave_idx_type iNumElemPerSlaveNode = 13;
     };

     template <>
     struct SelectElemPerSlaveNode<ShapeTria6H> {
          static constexpr octave_idx_type iNumElemPerSlaveNode = 13;
     };

     template <>
     struct SelectElemPerSlaveNode<ShapeTria10> {
          static constexpr octave_idx_type iNumElemPerSlaveNode = 13;
     };
};

class SurfToNodeConstrBase {
public:
     enum ConstraintType {
          CT_FIXED = 0,
          CT_SLIDING
     };

#if HAVE_NLOPT == 1
     enum ConstraintFlags: unsigned {
          CF_DEFAULT = 0x0u,
          CF_IGNORE_NODES_OUT_OF_RANGE = 0x1u,
          CF_ELEM_DOF_PRE_ALLOCATED = 0x2u
     };

     typedef ElemContact<std::complex<double>, ComplexMatrix> ElemContactType;
     typedef ElementBlock<ElemContactType> ElemBlockContactType;
     typedef ElementBlock<ElemJoint> ElemBlockJointType;
     typedef std::unique_ptr<ElemBlockJointType> ElemBlockJointPtr;
     typedef std::unique_ptr<ElemBlockContactType> ElemBlockContactPtr;

     static ConstraintType GetConstraintType(int constr) {
          switch (constr) {
          case CT_FIXED:
          case CT_SLIDING:
               return static_cast<ConstraintType>(constr);
          default:
               throw std::runtime_error("sfncon: invalid value for elements.sfncon{4|6|8}.constraint");
          }
     }

     static ConstraintType GetConstraintType(const Cell& ov, octave_idx_type j) {
          if (!ov.numel()) {
               return CT_FIXED;
          }

          const int constr = ov.xelem(j).int_value();

#if OCTAVE_MAJOR_VERSION < 6
          if (error_state) {
               throw std::runtime_error("sfncon: elements.sfncon{4|6|8}.constraint must be a scalar value");
          }
#endif

          return GetConstraintType(constr);
     }

     static octave_idx_type iGetNumDof(ConstraintType eType, DofMap::DomainType eDomain) {
          switch (eDomain) {
          case DofMap::DO_STRUCTURAL:
          case DofMap::DO_FLUID_STRUCT:
               return eType == CT_FIXED ? 3 : 1;
          case DofMap::DO_THERMAL:
          case DofMap::DO_ACOUSTICS:
               return 1;
          default:
               throw std::logic_error("sfncon: unknown value for dof_map.domain");
          }
     }

     static octave_idx_type iGetNumDof(const Cell& ov, octave_idx_type j, DofMap::DomainType eDomain) {
          return iGetNumDof(GetConstraintType(ov, j), eDomain);
     }

     static void BuildJoints(const Matrix& nodes,
                             const octave_scalar_map& elements,
                             const array<int32NDArray, DofMap::ELEM_TYPE_COUNT>& edof,
                             array<octave_idx_type, DofMap::ELEM_TYPE_COUNT>& dofelemid,
                             const ElementTypes::TypeInfo& oElemType,
                             const unsigned uConstraintFlags,
                             vector<std::unique_ptr<ElementBlockBase> >& rgElemBlocks,
                             DofMap::DomainType eDomain);

     static void BuildContacts(const Matrix& nodes,
                               const octave_scalar_map& elements,
                               const ElementTypes::TypeInfo& oElemType,
                               const unsigned uConstraintFlags,
                               vector<std::unique_ptr<ElementBlockBase> >& rgElemBlocks,
                               const DofMap::DomainType eDomain);
#else
     static const char szErrCompileWithNlopt[90];
#endif
};

#if HAVE_NLOPT == 0
const char SurfToNodeConstrBase::szErrCompileWithNlopt[] = "__mboct_fem_pkg__ must be compiled with nlopt in order to use element types sfncon{4|6|8}";
#endif

#if HAVE_NLOPT == 1
template <typename SHAPE_FUNC>
class SurfToNodeConstr: public SurfToNodeConstrBase {
public:
     static void
     BuildJoints(octave_idx_type& id,
                 const Matrix& X,
                 const int32NDArray& nidxmaster,
                 const int32NDArray& nidxslave,
                 const ColumnVector& maxdist,
                 const ConstraintType eType,
                 const unsigned uConstraintFlags,
                 const DofMap::DomainType eDomain,
                 vector<std::unique_ptr<ElementBlockBase> >& rgElemBlocks,
                 const double dScale) {

          ElementTypes::TypeId eElemType;

          switch (eDomain) {
          case DofMap::DO_STRUCTURAL:
          case DofMap::DO_FLUID_STRUCT:
               eElemType = ElementTypes::ELEM_JOINT;
               break;
          case DofMap::DO_THERMAL:
               eElemType = ElementTypes::ELEM_THERM_CONSTR;
               break;
          case DofMap::DO_ACOUSTICS:
               eElemType = ElementTypes::ELEM_ACOUSTIC_CONSTR;
               break;
          default:
               throw std::logic_error("sfncon: unsupported value for dof_map.domain");
          }

          const octave_idx_type iNumNodesSlave = nidxslave.numel();
          ElemBlockJointPtr pElemBlockPtr{new ElemBlockJointType{eElemType, iNumNodesSlave}};

          auto pfnInsertElemHook = [dScale, eType, eElemType, eDomain, &pElemBlockPtr]
               (const octave_idx_type id,
                const Matrix& X,
                const int32NDArray& nidxmaster,
                const int32NDArray& nidxslave,
                const unsigned uConstraintFlags,
                const octave_idx_type iSlaveNode,
                const octave_idx_type eidxmaster,
                const ColumnVector& rvopt,
                const ColumnVector& Xi) {
                                        InsertConstraint(id,
                                                         X,
                                                         nidxmaster,
                                                         nidxslave,
                                                         eType,
                                                         uConstraintFlags,
                                                         eDomain,
                                                         iSlaveNode,
                                                         pElemBlockPtr,
                                                         eElemType,
                                                         eidxmaster,
                                                         rvopt,
                                                         dScale);
                                   };

          FindSmallestDistance(id,
                               X,
                               nidxmaster,
                               nidxslave,
                               maxdist,
                               uConstraintFlags,
                               pfnInsertElemHook);

          rgElemBlocks.emplace_back(std::move(pElemBlockPtr));
     }

     static void
     BuildContacts(octave_idx_type& id,
                   const Matrix& X,
                   const int32NDArray& nidxmaster,
                   const int32NDArray& nidxslave,
                   const ColumnVector& maxdist,
                   const ConstraintType eType,
                   const unsigned uConstraintFlags,
                   const DofMap::DomainType eDomain,
                   vector<std::unique_ptr<ElementBlockBase>>& rgElemBlocks,
                   const octave_value& k) {

          const octave_idx_type iNumNodesSlave = nidxslave.numel();
          ElemBlockContactPtr pElemBlockPtr{new ElemBlockContactType{ElementTypes::ELEM_SPRING, iNumNodesSlave}};

          auto pfnInsertElemHook = [eType, eDomain, &k, &pElemBlockPtr]
               (const octave_idx_type id,
                const Matrix& X,
                const int32NDArray& nidxmaster,
                const int32NDArray& nidxslave,
                const unsigned uConstraintFlags,
                const octave_idx_type iSlaveNode,
                const octave_idx_type eidxmaster,
                const ColumnVector& rvopt,
                const ColumnVector& Xi) {
                                        InsertContact(id,
                                                      X,
                                                      nidxmaster,
                                                      nidxslave,
                                                      eType,
                                                      uConstraintFlags,
                                                      eDomain,
                                                      iSlaveNode,
                                                      pElemBlockPtr,
                                                      ElementTypes::ELEM_SPRING,
                                                      eidxmaster,
                                                      rvopt,
                                                      Xi,
                                                      k);
                                   };

          FindSmallestDistance(id,
                               X,
                               nidxmaster,
                               nidxslave,
                               maxdist,
                               uConstraintFlags,
                               pfnInsertElemHook);

          rgElemBlocks.emplace_back(std::move(pElemBlockPtr));
     }

private:
     typedef void InsertElemHook(const octave_idx_type id,
                                 const Matrix& X,
                                 const int32NDArray& nidxmaster,
                                 const int32NDArray& nidxslave,
                                 const unsigned uConstraintFlags,
                                 const octave_idx_type iSlaveNode,
                                 const octave_idx_type eidxmaster,
                                 const ColumnVector& rvopt,
                                 const ColumnVector& Xi);
     static void
     FindSmallestDistance(octave_idx_type& id,
                          const Matrix& X,
                          const int32NDArray& nidxmaster,
                          const int32NDArray& nidxslave,
                          const ColumnVector& maxdist,
                          const unsigned uConstraintFlags,
                          const std::function<InsertElemHook>& pfnInsertElemHook) {

          const octave_idx_type iNumNodesSlave = nidxslave.numel();
          const octave_idx_type iNumNodes = X.rows();
          const octave_idx_type iNumElem = nidxmaster.rows();

          FEM_ASSERT(X.columns() >= iNumDimNode);
          FEM_ASSERT(X.columns() == 6);
          FEM_ASSERT(nidxmaster.columns() == iNumNodesElem);
          FEM_ASSERT(nidxslave.rows() == maxdist.rows());

          ColumnVector Xs{iNumDimNode};

          vector<ElemIndexVector> eidxmaster(iNumNodesSlave);

          for (octave_idx_type k = 0; k < iNumNodesSlave; ++k) {
               for (octave_idx_type l = 0; l < iNumDimNode; ++l) {
                    Xs.xelem(l) = X.xelem(nidxslave.xelem(k).value() - 1 + iNumNodes * l);
               }

               for (octave_idx_type i = 0; i < nidxmaster.rows(); ++i) {
                    double dXmin = std::numeric_limits<double>::max();

                    for (octave_idx_type j = 0; j < iNumNodesElem; ++j) {
                         double dX = 0;

                         for (octave_idx_type l = 0; l < iNumDimNode; ++l) {
                              dX += std::pow(X.xelem(nidxmaster.xelem(i + iNumElem * j).value() - 1 + iNumNodes * l) - Xs.xelem(l), 2);
                         }

                         dXmin = std::min(dX, dXmin);
                    }

                    eidxmaster[k].insert(ElemIndexRecord(dXmin, i));
               }

               OCTAVE_QUIT;
          }

          ColumnVector Xm(iNumDimNode * iNumNodesElem);
          ColumnVector rv(iNumDir), rvopt(iNumDir, 0.);
          ColumnVector Xi(iNumDimNode), Xiopt(iNumDimNode, 0.);

          for (size_t i = 0; i < eidxmaster.size(); ++i) {
               for (octave_idx_type j = 0; j < iNumDimNode; ++j) {
                    Xs.xelem(j) = X.xelem(nidxslave.xelem(i).value() - 1, j);
               }

               double fopt = std::numeric_limits<double>::max();
               nlopt_result rcopt = NLOPT_FAILURE;
               octave_idx_type lopt = -1;
               rvopt.fill(0.);

               for (octave_idx_type l = 0; l < eidxmaster[i].size(); ++l) {
                    for (octave_idx_type j = 0; j < iNumNodesElem; ++j) {
                         for (octave_idx_type k = 0; k < iNumDimNode; ++k) {
                              Xm.xelem(j * iNumDimNode + k) = X.xelem(nidxmaster.xelem(eidxmaster[i][l].eidx, j).value() - 1, k);
                         }
                    }

                    rv.fill(0.);

                    double f = std::numeric_limits<double>::max();

                    nlopt_result rc = Solve(Xm, Xs, Xi, rv, f);

                    if (rc > 0 && f < fopt) {
                         rcopt = rc;
                         rvopt = rv;
                         fopt = f;
                         lopt = l;
                         Xiopt = Xi;
                    }
               }

               if (rcopt < 0 || fopt > maxdist.xelem(i) || lopt < 0) {
                    if (uConstraintFlags & CF_IGNORE_NODES_OUT_OF_RANGE) {
                         continue;
                    }

                    std::ostringstream os;

                    os << "sfncon: nlopt failed to project slave node #" << i + 1 << " (" << nidxslave(i).value()
                       << ") to element #";

                    octave_idx_type eidx = -1;

                    if (lopt >= 0) {
                         eidx = eidxmaster[i][lopt].eidx;

                         os << eidx + 1 << " (";

                         for (octave_idx_type k = 0; k < iNumNodesElem; ++k) {
                              os << nidxmaster.xelem(eidx, k).value() << ' ';
                         }

                         os << ')';
                    } else {
                         os << '?';
                    }

                    os << std::endl;
                    os << "status=" << rcopt
                       << " distance=" << fopt
                       << " maximum distance= " << maxdist(i)
                       << " position=" << rvopt.transpose()
                       << std::endl;

                    os << "Xs=";

                    for (octave_idx_type k = 0; k < iNumDimNode; ++k) {
                         os << X.xelem(nidxslave.xelem(i).value() - 1, k) << ' ';
                    }

                    os << std::endl;

                    os << "Xi=" << Xiopt.transpose() << std::endl;
                    os << "Xm=" << std::endl;

                    if (eidx >= 0) {
                         for (octave_idx_type k = 0; k < iNumNodesElem; ++k) {
                              for (octave_idx_type j = 0; j < iNumDimNode; ++j) {
                                   os << X.xelem(nidxmaster.xelem(eidx, k).value() - 1, j) << " ";
                              }
                              os << std::endl;
                         }
                    } else {
                         os << '?';
                    }

                    os << std::ends;

                    throw std::runtime_error(os.str());
               }

               FEM_TRACE("nlopt: i=" << i << " l=" << lopt << " rc=" << rcopt << " f=" << fopt << " maxdist=" << maxdist(i) << std::endl);

               pfnInsertElemHook(++id,
                                 X,
                                 nidxmaster,
                                 nidxslave,
                                 uConstraintFlags,
                                 i,
                                 eidxmaster[i][lopt].eidx,
                                 rvopt,
                                 Xi);
          }
     }

     static void InsertConstraint(const octave_idx_type id,
                                  const Matrix& X,
                                  const int32NDArray& nidxmaster,
                                  const int32NDArray& nidxslave,
                                  const ConstraintType eType,
                                  const unsigned uConstraintFlags,
                                  const DofMap::DomainType eDomain,
                                  const octave_idx_type iSlaveNode,
                                  const ElemBlockJointPtr& pElemBlock,
                                  const ElementTypes::TypeId eElemType,
                                  const octave_idx_type eidxmaster,
                                  const ColumnVector& rvopt,
                                  const double dScale) {
          octave_idx_type iNumDofNodeMax = -1;
          octave_idx_type iNumDofNodeConstr = -1;

          switch (eElemType) {
          case ElementTypes::ELEM_JOINT:
               iNumDofNodeMax = 6;
               iNumDofNodeConstr = 3;
               break;
          default:
               iNumDofNodeMax = 1;
               iNumDofNodeConstr = 1;
          }

          const octave_idx_type iNumNodes = X.rows();

          Matrix Hf(iNumDofNodeConstr, iNumDofNodeConstr * iNumNodesElem);

          switch (eElemType) {
          case ElementTypes::ELEM_JOINT:
               SHAPE_FUNC::VectorInterpMatrix(rvopt, Hf);
               break;
          default:
               SHAPE_FUNC::ScalarInterpMatrix(rvopt, Hf, 0);
          }

          Matrix C(iNumDofNodeConstr, iNumDofNodeMax * (iNumNodesElem + 1), 0.);

          for (octave_idx_type j = 0; j < iNumDofNodeConstr; ++j) {
               for (octave_idx_type k = 0; k < iNumDofNodeConstr; ++k) {
                    C.xelem(k + j * iNumDofNodeConstr) = (k == j) ? 1. : 0.;
               }
          }

          for (octave_idx_type j = 0; j < iNumNodesElem; ++j) {
               for (octave_idx_type k = 0; k < iNumDofNodeConstr; ++k) {
                    for (octave_idx_type l = 0; l < iNumDofNodeConstr; ++l) {
                         C.xelem(l + iNumDofNodeConstr * ((j + 1) * iNumDofNodeMax + k)) = -Hf.xelem(l + iNumDofNodeConstr * (j * iNumDofNodeConstr + k));
                    }
               }
          }

          int32NDArray enodes(dim_vector(iNumNodesElem + 1, 1));

          enodes.xelem(0) = nidxslave.xelem(iSlaveNode).value();

          for (octave_idx_type j = 0; j < iNumNodesElem; ++j) {
               enodes.xelem(j + 1) = nidxmaster.xelem(eidxmaster, j).value();
          }

          const octave_idx_type iNumCoord = X.columns();

          FEM_ASSERT(iNumCoord == 6);

          Matrix Xe(iNumCoord, iNumNodesElem + 1);

          for (octave_idx_type j = 0; j < iNumCoord; ++j) {
               Xe.xelem(j) = X.xelem(nidxslave.xelem(iSlaveNode).value() - 1, j);
          }

          for (octave_idx_type j = 0; j < iNumNodesElem; ++j) {
               for (octave_idx_type k = 0; k < iNumCoord; ++k) {
                    Xe.xelem(k + iNumCoord * (j + 1)) = X.xelem(nidxmaster.xelem(eidxmaster, j).value() - 1 + iNumNodes * k);
               }
          }

          if (eType == CT_SLIDING) {
               if (eElemType != ElementTypes::ELEM_JOINT) {
                    throw std::logic_error("sfncon: unsupported value for dof_map.domain");
               }

               Matrix dHf_dr(iNumDimNode, iNumDofNodeConstr * iNumNodesElem);
               Matrix dHf_ds(iNumDimNode, iNumDofNodeConstr * iNumNodesElem);
               ColumnVector n1(iNumDimNode), n2(iNumDimNode);
               ColumnVector n(iNumDimNode);
               RowVector nC(C.columns(), 0.);

               SHAPE_FUNC::VectorInterpMatrixDerR(rvopt, dHf_dr);
               SHAPE_FUNC::VectorInterpMatrixDerS(rvopt, dHf_ds);

               SurfaceTangentVector(Xe, dHf_dr, n1);
               SurfaceTangentVector(Xe, dHf_ds, n2);
               SurfaceElement::SurfaceNormalVectorUnit(n1, n2, n);

               for (octave_idx_type j = 0; j < C.columns(); ++j) {
                    for (octave_idx_type k = 0; k < iNumDofNodeConstr; ++k) {
                         nC.xelem(j) += n.xelem(k) * C.xelem(k + iNumDofNodeConstr * j);
                    }
               }

               FEM_TRACE("n=" << n << std::endl);
               FEM_TRACE("C=" << C << std::endl);
               FEM_TRACE("n.'*C=" << nC << std::endl);

               pElemBlock->Insert(id, Xe, nullptr, enodes, nC, Matrix{1, 0}, eDomain, dScale);
          } else {
               pElemBlock->Insert(id, Xe, nullptr, enodes, C, Matrix{iNumDofNodeConstr, 0}, eDomain, dScale);
          }
     }

     static void InsertContact(const octave_idx_type id,
                               const Matrix& X,
                               const int32NDArray& nidxmaster,
                               const int32NDArray& nidxslave,
                               const ConstraintType eType,
                               const unsigned uConstraintFlags,
                               const DofMap::DomainType eDomain,
                               const octave_idx_type iSlaveNode,
                               const ElemBlockContactPtr& pElemBlock,
                               const ElementTypes::TypeId eElemType,
                               const octave_idx_type eidxmaster,
                               const ColumnVector& rvopt,
                               const ColumnVector& Xi,
                               const octave_value& k) {
          constexpr octave_idx_type iNumDofNodeConstr = 3;
          const octave_idx_type iNumNodes = X.rows();

          Matrix Hf(iNumDofNodeConstr, iNumDofNodeConstr * iNumNodesElem);

          SHAPE_FUNC::VectorInterpMatrix(rvopt, Hf);

          int32NDArray enodes(dim_vector(iNumNodesElem + 1, 1));

          enodes.xelem(0) = nidxslave.xelem(iSlaveNode).value();

          for (octave_idx_type j = 0; j < iNumNodesElem; ++j) {
               enodes.xelem(j + 1) = nidxmaster.xelem(eidxmaster, j).value();
          }

          const octave_idx_type iNumCoord = X.columns();

          FEM_ASSERT(iNumCoord == 6);

          Matrix Xe(iNumCoord, iNumNodesElem + 1);

          for (octave_idx_type j = 0; j < iNumCoord; ++j) {
               Xe.xelem(j) = X.xelem(nidxslave.xelem(iSlaveNode).value() - 1, j);
          }

          for (octave_idx_type j = 0; j < iNumNodesElem; ++j) {
               for (octave_idx_type k = 0; k < iNumCoord; ++k) {
                    Xe.xelem(k + iNumCoord * (j + 1)) = X.xelem(nidxmaster.xelem(eidxmaster, j).value() - 1 + iNumNodes * k);
               }
          }

          Matrix dHf_dr(iNumDimNode, iNumDofNodeConstr * iNumNodesElem);
          Matrix dHf_ds(iNumDimNode, iNumDofNodeConstr * iNumNodesElem);
          ColumnVector n1(iNumDimNode), n2(iNumDimNode);
          ColumnVector n3(iNumDimNode);

          SHAPE_FUNC::VectorInterpMatrixDerR(rvopt, dHf_dr);
          SHAPE_FUNC::VectorInterpMatrixDerS(rvopt, dHf_ds);

          SurfaceTangentVector(Xe, dHf_dr, n1);
          SurfaceTangentVector(Xe, dHf_ds, n2);
          SurfaceElement::SurfaceNormalVector(n1, n2, n3);
          const double dA = SurfaceElement::JacobianDet(n3);

          Matrix R(iNumDimNode, iNumDimNode);

          for (octave_idx_type i = 0; i < iNumDimNode; ++i) {
               R.xelem(i, 0) = n1.xelem(i);
               R.xelem(i, 1) = n2.xelem(i);
               R.xelem(i, 2) = n3.xelem(i);
          }

          for (octave_idx_type j = 0; j < iNumDimNode; ++j) {
               double n = 0.;

               for (octave_idx_type i = 0; i < iNumDimNode; ++i) {
                    n += std::pow(R.xelem(i, j), 2);
               }

               n = sqrt(n);

               for (octave_idx_type i = 0; i < iNumDimNode; ++i) {
                    R.xelem(i, j) /= n;
               }
          }

          double h0 = 0.;

          for (octave_idx_type i = 0; i < iNumDimNode; ++i) {
               h0 += R.xelem(i, 2) * (Xe.xelem(i) - Xi.xelem(i));
          }

          pElemBlock->Insert(id, Xe, nullptr, enodes, R, Hf, k, h0, dA);
     }

     static void SurfaceTangentVector(const Matrix& X, const Matrix& dHf, ColumnVector& n) {
          FEM_ASSERT(dHf.rows() == 3);

          const octave_idx_type iNumCoord = X.rows();
          const octave_idx_type iNumNodes = X.columns();

          for (octave_idx_type i = 0; i < 3; ++i) {
               double ni = 0.;

               for (octave_idx_type j = 0; j < iNumNodes - 1; ++j) {
                    for (octave_idx_type k = 0; k < 3; ++k) {
                         ni += dHf.xelem(i + 3 * (j * 3 + k)) * X.xelem(k + iNumCoord * (j + 1));
                    }
               }

               n.xelem(i) = ni;
          }
     }

     static nlopt_result
     Solve(const ColumnVector& Xm,
           const ColumnVector& Xs,
           ColumnVector& Xi,
           ColumnVector& rv,
           double& f) {

          constexpr double dTolX = std::pow(std::numeric_limits<double>::epsilon(), 0.8);

          FuncData oFuncData(Xm, Xs);

          nlopt_set_min_objective(oFuncData.opt, &SurfToNodeConstr::Objective, &oFuncData);

          if (SHAPE_FUNC::iGetNumEqualityConstr()) {
              if (nlopt_add_equality_constraint(oFuncData.opt, &SurfToNodeConstr::EqualityConstr, &oFuncData, dTolX) < 0) {
                  throw std::runtime_error("sfncon: nlopt_add_equality_constraint failed");
              }
          }

          if (SHAPE_FUNC::iGetNumInequalityConstr()) {
              if (nlopt_add_inequality_constraint(oFuncData.opt, &SurfToNodeConstr::InequalityConstr, &oFuncData, dTolX) < 0) {
                  throw std::runtime_error("sfncon: nlopt_add_inequality_constraint failed");
              }
          }

          Matrix dX(iNumDimNode, 2);

          for (octave_idx_type i = 0; i < iNumDimNode; ++i) {
               dX.xelem(i) = std::numeric_limits<double>::max();
               dX.xelem(i + iNumDimNode) = -dX.xelem(i);

               for (octave_idx_type j = 0; j < iNumNodesElem; ++j) {
                    dX.xelem(i) = std::min(dX.xelem(i), Xm.xelem(j * iNumDimNode + i));
                    dX.xelem(i + iNumDimNode) = std::max(dX.xelem(i + iNumDimNode), Xm.xelem(j * iNumDimNode + i));
               }
          }

          double dTolF = 0;

          for (octave_idx_type i = 0; i < iNumDimNode; ++i) {
               dTolF = std::max(dTolF, dTolX * (dX.xelem(i + iNumDimNode) - dX.xelem(i)));
          }

          dTolF *= dTolF;

          ColumnVector rmin(iNumDir);
          ColumnVector rmax(iNumDir);

          SHAPE_FUNC::GetElemLimits(rmin, rmax);

          nlopt_set_lower_bounds(oFuncData.opt, rmin.fortran_vec());
          nlopt_set_upper_bounds(oFuncData.opt, rmax.fortran_vec());
          nlopt_set_maxeval(oFuncData.opt, 10000);
          nlopt_set_xtol_abs1(oFuncData.opt, dTolX);
          nlopt_set_ftol_abs(oFuncData.opt, dTolF);
          auto rc = nlopt_optimize(oFuncData.opt, rv.fortran_vec(), &f);
          oFuncData.oSNCO.Position(rv, Xi, oFuncData.Hf);

          f = sqrt(f);

          return rc;
     }

     static constexpr octave_idx_type iNumDimNode = SHAPE_FUNC::iGetNumDofNode();
     static constexpr octave_idx_type iNumDir = SHAPE_FUNC::iGetNumDirections();
     static constexpr octave_idx_type iNumNodesElem = SHAPE_FUNC::iGetNumNodes();

     SurfToNodeConstr(const ColumnVector& Xm, const ColumnVector& Xs)
          :Xm(Xm),
           Xs(Xs) {

          FEM_ASSERT(Xm.rows() == iNumNodesElem * Xs.rows());
          FEM_ASSERT(Xs.rows() == iNumDimNode);
     }

     static double Objective(unsigned n, const double x[], double gradient[], void* pData) {
          auto pFuncData = static_cast<FuncData*>(pData);

          Matrix& Hf = pFuncData->Hf;
          ColumnVector& rv = pFuncData->rv;
          ColumnVector& F = pFuncData->f;
          const SurfToNodeConstr& oSNCO = pFuncData->oSNCO;
          const octave_idx_type N = rv.rows();

          FEM_ASSERT(n == N);

          for (octave_idx_type i = 0; i < N; ++i) {
               rv.xelem(i) = x[i];
          }

          const double f0 = oSNCO.Objective(rv, F, Hf);

          if (gradient) {
               constexpr double delta = std::pow(std::numeric_limits<double>::epsilon(), 0.5);
               constexpr octave_idx_type M = 4;
               static constexpr double dx[M] = {-2. * delta, -delta, delta, 2. * delta};
               std::array<double, M> fi;

               for (octave_idx_type i = 0; i < N; ++i) {
                    for (octave_idx_type j = 0; j < M; ++j) {
                         rv.xelem(i) = x[i] + dx[j];
                         fi[j] = oSNCO.Objective(rv, F, Hf);
                    }

                    gradient[i] = (fi[0] - 8. * fi[1] + 8. * fi[2] - fi[3]) / (12. * delta);
                    rv.xelem(i) = x[i];
               }
          }

          return f0;
     }

     static double EqualityConstr(unsigned n, const double x[], double gradient[], void* pData) {
          auto pFuncData = static_cast<FuncData*>(pData);

          ColumnVector& rv = pFuncData->rv;
          const SurfToNodeConstr& oSNCO = pFuncData->oSNCO;
          const octave_idx_type N = rv.rows();

          FEM_ASSERT(n == N);

          for (octave_idx_type i = 0; i < N; ++i) {
               rv.xelem(i) = x[i];
          }

          const double f0 = oSNCO.EqualityConstr(rv);

          if (gradient) {
               constexpr double delta = std::pow(std::numeric_limits<double>::epsilon(), 0.5);
               constexpr octave_idx_type M = 4;
               static constexpr double dx[M] = {-2. * delta, -delta, delta, 2. * delta};
               std::array<double, M> fi;

               for (octave_idx_type i = 0; i < N; ++i) {
                    for (octave_idx_type j = 0; j < M; ++j) {
                         rv.xelem(i) = x[i] + dx[j];
                         fi[j] = oSNCO.EqualityConstr(rv);
                    }

                    gradient[i] = (fi[0] - 8. * fi[1] + 8. * fi[2] - fi[3]) / (12. * delta);
                    rv.xelem(i) = x[i];
               }
          }

          return f0;
     }

     static double InequalityConstr(unsigned n, const double x[], double gradient[], void* pData) {
          auto pFuncData = static_cast<FuncData*>(pData);

          ColumnVector& rv = pFuncData->rv;
          const SurfToNodeConstr& oSNCO = pFuncData->oSNCO;
          const octave_idx_type N = rv.rows();

          FEM_ASSERT(n == N);

          for (octave_idx_type i = 0; i < N; ++i) {
               rv.xelem(i) = x[i];
          }

          const double f0 = oSNCO.InequalityConstr(rv);

          if (gradient) {
               constexpr double delta = std::pow(std::numeric_limits<double>::epsilon(), 0.5);
               constexpr octave_idx_type M = 4;
               static constexpr double dx[M] = {-2. * delta, -delta, delta, 2. * delta};
               std::array<double, M> fi;

               for (octave_idx_type i = 0; i < N; ++i) {
                    for (octave_idx_type j = 0; j < M; ++j) {
                         rv.xelem(i) = x[i] + dx[j];
                         fi[j] = oSNCO.InequalityConstr(rv);
                    }

                    gradient[i] = (fi[0] - 8. * fi[1] + 8. * fi[2] - fi[3]) / (12. * delta);
                    rv.xelem(i) = x[i];
               }
          }

          return f0;
     }

     double Objective(const ColumnVector& rv, ColumnVector& f, Matrix& Hf) const {
          FEM_ASSERT(f.rows() == iNumDimNode);
          FEM_ASSERT(f.rows() == Xs.rows());

          Position(rv, f, Hf);

          for (octave_idx_type i = 0; i < f.rows(); ++i) {
               f.xelem(i) -= Xs.xelem(i);
          }

          double ftot = 0.;

          for (octave_idx_type i = 0; i < f.rows(); ++i) {
               ftot += f.xelem(i) * f.xelem(i);
          }

          OCTAVE_QUIT;
#if DEBUG > 2
          FEM_TRACE("f=" << sqrt(ftot) << std::endl);
#endif
          return ftot;
     }

     static double EqualityConstr(const ColumnVector& rv) {
          FEM_ASSERT(SHAPE_FUNC::iGetNumEqualityConstr());

          double f = SHAPE_FUNC::EqualityConstr(rv);

          OCTAVE_QUIT;

          return f;
     }

     static double InequalityConstr(const ColumnVector& rv) {
          FEM_ASSERT(SHAPE_FUNC::iGetNumInequalityConstr());

          double f = SHAPE_FUNC::InequalityConstr(rv);

          OCTAVE_QUIT;

          return f;
     }

     void Position(const ColumnVector& rv, ColumnVector& Xi, Matrix& Hf) const {
          FEM_ASSERT(rv.rows() == iNumDir);
          FEM_ASSERT(Xi.rows() == iNumDimNode);
          FEM_ASSERT(Hf.rows() == iNumDimNode);
          FEM_ASSERT(Hf.columns() == Xi.rows() * iNumNodesElem);
          FEM_ASSERT(Xm.rows() == Hf.columns());

          SHAPE_FUNC::VectorInterpMatrix(rv, Hf);

          for (octave_idx_type i = 0; i < Xi.rows(); ++i) {
               Xi.xelem(i) = 0;
          }

          for (octave_idx_type j = 0; j < Hf.columns(); ++j) {
               for (octave_idx_type i = 0; i < iNumDimNode; ++i) {
                    Xi.xelem(i) += Hf.xelem(i + iNumDimNode * j) * Xm.xelem(j);
               }
          }
     }

private:
     const ColumnVector Xm, Xs;

     struct FuncData {
          FuncData(const ColumnVector& Xm, const ColumnVector& Xs)
               :oSNCO(Xm, Xs),
                opt(nlopt_create(NLOPT_LD_SLSQP, iNumDir)),
                Hf(iNumDimNode, iNumNodesElem * iNumDimNode),
                rv(iNumDir, 0.),
                f(iNumDimNode) {
          }

          ~FuncData() {
               nlopt_destroy(opt);
          }

          const SurfToNodeConstr oSNCO;
          nlopt_opt opt;
          Matrix Hf;
          ColumnVector rv, f;
     };

     static constexpr octave_idx_type iNumElemPerSlaveNode = shape_func_util::SelectElemPerSlaveNode<SHAPE_FUNC>::iNumElemPerSlaveNode;

     struct ElemIndexRecord {
          ElemIndexRecord(double dX = std::numeric_limits<double>::max(), octave_idx_type eidx = -1)
               :dX(dX),
                eidx(eidx) {
          }

          bool operator<(const ElemIndexRecord& rhs) const {
               return dX < rhs.dX;
          }

          double dX;
          octave_idx_type eidx;
     };

     class ElemIndexVector {
     public:
          ElemIndexVector()
               :isize(0) {
          }

          octave_idx_type size() const {
               return isize;
          }

          void insert(const ElemIndexRecord& oRec) {
               decltype(data.begin()) newrecord;

               if (size() >= iNumElemPerSlaveNode) {
                    newrecord = std::max_element(data.begin(), data.end());
               } else {
                    newrecord = data.begin() + size();
                    ++isize;
               }

               if (oRec.dX <= newrecord->dX) {
                    *newrecord = oRec;
               }
          }

          ElemIndexRecord& operator[](octave_idx_type i) {
               FEM_ASSERT(i >= 0);
               FEM_ASSERT(i < size());
               return data[i];
          }

     private:
          array<ElemIndexRecord, iNumElemPerSlaveNode> data;
          octave_idx_type isize;
     };
};

void SurfToNodeConstrBase::BuildJoints(const Matrix& nodes,
                                       const octave_scalar_map& elements,
                                       const array<int32NDArray, DofMap::ELEM_TYPE_COUNT>& edof,
                                       array<octave_idx_type, DofMap::ELEM_TYPE_COUNT>& dofelemid,
                                       const ElementTypes::TypeInfo& oElemType,
                                       const unsigned uConstraintFlags,
                                       vector<std::unique_ptr<ElementBlockBase> >& rgElemBlocks,
                                       const DofMap::DomainType eDomain) {
     const auto iter_elem = elements.seek(oElemType.name);

     if (iter_elem == elements.end()) {
          return;
     }

     const octave_map s_elem(elements.contents(iter_elem).map_value());

#if OCTAVE_MAJOR_VERSION < 6
     if (error_state) {
          throw std::runtime_error("sfncon: elements.sfncon{4|6|8} must be a struct array");
     }
#endif

     const auto iter_nidxmaster = s_elem.seek("master");

     if (iter_nidxmaster == s_elem.end()) {
          throw std::runtime_error("sfncon: elements.sfncon{4|6|8}.master not defined");
     }

     const Cell ov_nidxmaster = s_elem.contents(iter_nidxmaster);

     const auto iter_nidxslave = s_elem.seek("slave");

     if (iter_nidxslave == s_elem.end()) {
          throw std::runtime_error("sfncon: elements.sfncon{4|6|8}.slave not defined");
     }

     const Cell ov_nidxslave = s_elem.contents(iter_nidxslave);

     const auto iter_maxdist = s_elem.seek("maxdist");

     if (iter_maxdist == s_elem.end()) {
          throw std::runtime_error("sfncon: elements.sfncon{4|6|8}.maxdist not defined");
     }

     const Cell ov_maxdist = s_elem.contents(iter_maxdist);

     const auto iter_constr = s_elem.seek("constraint");

     Cell ov_constr;

     if (iter_constr != s_elem.end()) {
          ov_constr = s_elem.contents(iter_constr);
     }

     Cell ov_Scale;

     const auto iter_Scale = s_elem.seek("scale");

     if (iter_Scale != s_elem.end()) {
          ov_Scale = s_elem.contents(iter_Scale);
     }

     FEM_ASSERT(ov_nidxslave.numel() == s_elem.numel());
     FEM_ASSERT(ov_nidxmaster.numel() == s_elem.numel());
     FEM_ASSERT(ov_maxdist.numel() == s_elem.numel());
     FEM_ASSERT(ov_constr.numel() == 0 || ov_constr.numel() == s_elem.numel());
     FEM_ASSERT(ov_Scale.numel() == 0 || ov_Scale.numel() == s_elem.numel());

     for (octave_idx_type l = 0; l < s_elem.numel(); ++l) {
          const int32NDArray nidxmaster = ov_nidxmaster(l).int32_array_value();

#if OCTAVE_MAJOR_VERSION < 6
          if (error_state) {
               throw std::runtime_error("sfncon: elements.sfncon{4|6|8}.master must be an integer array");
          }
#endif

          const int32NDArray nidxslave = ov_nidxslave(l).int32_array_value();

#if OCTAVE_MAJOR_VERSION < 6
          if (error_state) {
               throw std::runtime_error("sfncon: elements.sfncon{4|6|8}.slave must be an integer array");
          }
#endif

          ColumnVector maxdist = ov_maxdist(l).column_vector_value();

#if OCTAVE_MAJOR_VERSION < 6
          if (error_state) {
               throw std::runtime_error("sfncon: elements.sfncon{4|6|8}.maxdist must be a column vector");
          }
#endif

          const auto eConstrType = SurfToNodeConstrBase::GetConstraintType(ov_constr, l);

          octave_idx_type iNumNodesElem = -1;

          switch (oElemType.type) {
          case ElementTypes::ELEM_SFNCON4:
               iNumNodesElem = ShapeIso4::iGetNumNodes();
               break;
          case ElementTypes::ELEM_SFNCON6:
               iNumNodesElem = ShapeTria6::iGetNumNodes();
               break;
          case ElementTypes::ELEM_SFNCON6H:
               iNumNodesElem = ShapeTria6H::iGetNumNodes();
               break;
          case ElementTypes::ELEM_SFNCON8:
               iNumNodesElem = ShapeQuad8::iGetNumNodes();
               break;
          case ElementTypes::ELEM_SFNCON8R:
               iNumNodesElem = ShapeQuad8r::iGetNumNodes();
               break;
          case ElementTypes::ELEM_SFNCON9:
               iNumNodesElem = ShapeQuad9::iGetNumNodes();
               break;
          case ElementTypes::ELEM_SFNCON10:
               iNumNodesElem = ShapeTria10::iGetNumNodes();
               break;
          default:
               FEM_ASSERT(false);
          }

          if (nidxmaster.ndims() != 2 || nidxmaster.rows() < 1 || nidxmaster.columns() != iNumNodesElem) {
               throw std::runtime_error("sfncon: elements.sfncon{4|6|8}.master must be an nx{4|6|8} array");
          }

          if (nidxslave.ndims() != 2 || nidxslave.rows() < 1 || nidxslave.columns() != 1) {
               throw std::runtime_error("snfcon: elements.sfncon{4|6|8}.slave must be an nx1 array");
          }

          if (maxdist.rows() == 1 && nidxslave.rows() > 1) {
               const double maxdistval = maxdist(0);
               maxdist.resize(nidxslave.rows(), maxdistval);
          }

          if (maxdist.rows() != nidxslave.rows()) {
               throw std::runtime_error("sfncon: elements.sfncon{4|6|8}.maxdist must have "
                                        "the same dimensions like elements.sfncon{4|6|8}.slave");
          }

          for (octave_idx_type i = 0; i < nidxmaster.rows(); ++i) {
               for (octave_idx_type j = 0; j < nidxmaster.columns(); ++j) {
                    if (nidxmaster(i, j).value() < 1 || nidxmaster(i, j).value() > nodes.rows()) {
                         throw std::runtime_error("sfncon: elements.sfncon{4|6|8}.master: node index out of range");
                    }
               }
          }

          for (octave_idx_type i = 0; i < nidxslave.rows(); ++i) {
               if (nidxslave(i).value() < 1 || nidxslave(i).value() > nodes.rows()) {
                    throw std::runtime_error("sfncon: elements.sfncon{4|6|8}.slave: node index out of range");
               }
          }

          const double dScale = ov_Scale.numel() ? ov_Scale.xelem(l).scalar_value() : 1.;

#if OCTAVE_MAJOR_VERSION < 6
          if (error_state) {
               throw std::runtime_error("sfncon: mesh.elements.sfncon{4|6|8}.scale must be a real scalar");
          }
#endif

          switch (oElemType.type) {
          case ElementTypes::ELEM_SFNCON4:
               SurfToNodeConstr<ShapeIso4>::BuildJoints(dofelemid[oElemType.dof_type],
                                                        nodes,
                                                        nidxmaster,
                                                        nidxslave,
                                                        maxdist,
                                                        eConstrType,
                                                        uConstraintFlags,
                                                        eDomain,
                                                        rgElemBlocks,
                                                        dScale);
               break;
          case ElementTypes::ELEM_SFNCON6:
               SurfToNodeConstr<ShapeTria6>::BuildJoints(dofelemid[oElemType.dof_type],
                                                         nodes,
                                                         nidxmaster,
                                                         nidxslave,
                                                         maxdist,
                                                         eConstrType,
                                                         uConstraintFlags,
                                                         eDomain,
                                                         rgElemBlocks,
                                                         dScale);
               break;
          case ElementTypes::ELEM_SFNCON6H:
               SurfToNodeConstr<ShapeTria6H>::BuildJoints(dofelemid[oElemType.dof_type],
                                                          nodes,
                                                          nidxmaster,
                                                          nidxslave,
                                                          maxdist,
                                                          eConstrType,
                                                          uConstraintFlags,
                                                          eDomain,
                                                          rgElemBlocks,
                                                          dScale);
               break;
          case ElementTypes::ELEM_SFNCON8:
               SurfToNodeConstr<ShapeQuad8>::BuildJoints(dofelemid[oElemType.dof_type],
                                                         nodes,
                                                         nidxmaster,
                                                         nidxslave,
                                                         maxdist,
                                                         eConstrType,
                                                         uConstraintFlags,
                                                         eDomain,
                                                         rgElemBlocks,
                                                         dScale);
               break;
          case ElementTypes::ELEM_SFNCON8R:
               SurfToNodeConstr<ShapeQuad8r>::BuildJoints(dofelemid[oElemType.dof_type],
                                                          nodes,
                                                          nidxmaster,
                                                          nidxslave,
                                                          maxdist,
                                                          eConstrType,
                                                          uConstraintFlags,
                                                          eDomain,
                                                          rgElemBlocks,
                                                          dScale);
               break;
          case ElementTypes::ELEM_SFNCON9:
               SurfToNodeConstr<ShapeQuad9>::BuildJoints(dofelemid[oElemType.dof_type],
                                                         nodes,
                                                         nidxmaster,
                                                         nidxslave,
                                                         maxdist,
                                                         eConstrType,
                                                         uConstraintFlags,
                                                         eDomain,
                                                         rgElemBlocks,
                                                         dScale);
               break;
          case ElementTypes::ELEM_SFNCON10:
               SurfToNodeConstr<ShapeTria10>::BuildJoints(dofelemid[oElemType.dof_type],
                                                          nodes,
                                                          nidxmaster,
                                                          nidxslave,
                                                          maxdist,
                                                          eConstrType,
                                                          uConstraintFlags,
                                                          eDomain,
                                                          rgElemBlocks,
                                                          dScale);
               break;
          default:
               FEM_ASSERT(false);
          }

          if ((uConstraintFlags & CF_ELEM_DOF_PRE_ALLOCATED) && dofelemid[oElemType.dof_type] > edof[oElemType.dof_type].rows()) {
               throw std::runtime_error("sfncon: dof_map.edof is not consistent with elements");
          }
     }
}

void SurfToNodeConstrBase::BuildContacts(const Matrix& nodes,
                                         const octave_scalar_map& elements,
                                         const ElementTypes::TypeInfo& oElemType,
                                         const unsigned uConstraintFlags,
                                         vector<std::unique_ptr<ElementBlockBase> >& rgElemBlocks,
                                         const DofMap::DomainType eDomain) {
     const auto iter_elem = elements.seek(oElemType.name);

     if (iter_elem == elements.end()) {
          return;
     }

     const auto eConstrType = CT_FIXED;

     const octave_map s_elem(elements.contents(iter_elem).map_value());

#if OCTAVE_MAJOR_VERSION < 6
     if (error_state) {
          throw std::runtime_error("sfncon: elements.sfncon{4|6|8}s must be a struct array");
     }
#endif

     const auto iter_nidxmaster = s_elem.seek("master");

     if (iter_nidxmaster == s_elem.end()) {
          throw std::runtime_error("sfncon: elements.sfncon{4|6|8}s.master not defined");
     }

     const Cell ov_nidxmaster = s_elem.contents(iter_nidxmaster);

     const auto iter_nidxslave = s_elem.seek("slave");

     if (iter_nidxslave == s_elem.end()) {
          throw std::runtime_error("sfncon: elements.sfncon{4|6|8}s.slave not defined");
     }

     const Cell ov_nidxslave = s_elem.contents(iter_nidxslave);

     const auto iter_maxdist = s_elem.seek("maxdist");

     if (iter_maxdist == s_elem.end()) {
          throw std::runtime_error("sfncon: elements.sfncon{4|6|8}s.maxdist not defined");
     }

     const Cell ov_maxdist = s_elem.contents(iter_maxdist);

     const auto iter_k = s_elem.seek("k");

     if (iter_k == s_elem.end()) {
          throw std::runtime_error("sfncon: elements.sfncon{4|6|8}s.k not defined");
     }

     Cell ov_k = s_elem.contents(iter_k);

     // const auto iter_sigma0 = s_elem.seek("sigma0");

     // if (iter_sigma0 == s_elem.end()) {
     //      throw std::runtime_error("sfncon: elements.sfncon{4|6|8}s.sigma0 not defined");
     // }

     // Cell ov_sigma0 = s_elem.contents(iter_sigma0);

     // const auto iter_sigma_delta = s_elem.seek("sigma_delta");

     // if (iter_sigma_delta == s_elem.end()) {
     //      throw std::runtime_error("sfncon: elements.sfncon{4|6|8}s.sigma_delta not defined");
     // }

     // Cell ov_sigma_delta = s_elem.contents(iter_sigma_delta);

     FEM_ASSERT(ov_nidxslave.numel() == s_elem.numel());
     FEM_ASSERT(ov_nidxmaster.numel() == s_elem.numel());
     FEM_ASSERT(ov_maxdist.numel() == s_elem.numel());

     octave_idx_type id = 0;

     for (octave_idx_type l = 0; l < s_elem.numel(); ++l) {
          const int32NDArray nidxmaster = ov_nidxmaster(l).int32_array_value();

#if OCTAVE_MAJOR_VERSION < 6
          if (error_state) {
               throw std::runtime_error("sfncon: elements.sfncon{4|6|8}s.master must be an integer array");
          }
#endif

          const int32NDArray nidxslave = ov_nidxslave(l).int32_array_value();

#if OCTAVE_MAJOR_VERSION < 6
          if (error_state) {
               throw std::runtime_error("sfncon: elements.sfncon{4|6|8}.slave must be an integer array");
          }
#endif

          ColumnVector maxdist = ov_maxdist(l).column_vector_value();

#if OCTAVE_MAJOR_VERSION < 6
          if (error_state) {
               throw std::runtime_error("sfncon: elements.sfncon{4|6|8}.maxdist must be a column vector");
          }
#endif
          octave_idx_type iNumNodesElem = -1;

          switch (oElemType.type) {
          case ElementTypes::ELEM_SFNCON4S:
               iNumNodesElem = ShapeIso4::iGetNumNodes();
               break;
          case ElementTypes::ELEM_SFNCON6S:
               iNumNodesElem = ShapeTria6::iGetNumNodes();
               break;
          case ElementTypes::ELEM_SFNCON6HS:
               iNumNodesElem = ShapeTria6H::iGetNumNodes();
               break;
          case ElementTypes::ELEM_SFNCON8S:
               iNumNodesElem = ShapeQuad8::iGetNumNodes();
               break;
          case ElementTypes::ELEM_SFNCON8RS:
               iNumNodesElem = ShapeQuad8r::iGetNumNodes();
               break;
          case ElementTypes::ELEM_SFNCON9S:
               iNumNodesElem = ShapeQuad9::iGetNumNodes();
               break;
          case ElementTypes::ELEM_SFNCON10S:
               iNumNodesElem = ShapeTria10::iGetNumNodes();
               break;
          default:
               FEM_ASSERT(false);
          }

          if (nidxmaster.ndims() != 2 || nidxmaster.rows() < 1 || nidxmaster.columns() != iNumNodesElem) {
               throw std::runtime_error("sfncon: elements.sfncon{4|6|8}s.master must be an nx{4|6|8} array");
          }

          if (nidxslave.ndims() != 2 || nidxslave.rows() < 1 || nidxslave.columns() != 1) {
               throw std::runtime_error("snfcon: elements.sfncon{4|6|8}s.slave must be an nx1 array");
          }

          if (maxdist.rows() == 1 && nidxslave.rows() > 1) {
               const double maxdistval = maxdist(0);
               maxdist.resize(nidxslave.rows(), maxdistval);
          }

          if (maxdist.rows() != nidxslave.rows()) {
               throw std::runtime_error("sfncon: elements.sfncon{4|6|8}s.maxdist must have "
                                        "the same dimensions like elements.sfncon{4|6|8}s.slave");
          }

          for (octave_idx_type i = 0; i < nidxmaster.rows(); ++i) {
               for (octave_idx_type j = 0; j < nidxmaster.columns(); ++j) {
                    if (nidxmaster(i, j).value() < 1 || nidxmaster(i, j).value() > nodes.rows()) {
                         throw std::runtime_error("sfncon: elements.sfncon{4|6|8}s.master: node index out of range");
                    }
               }
          }

          for (octave_idx_type i = 0; i < nidxslave.rows(); ++i) {
               if (nidxslave(i).value() < 1 || nidxslave(i).value() > nodes.rows()) {
                    throw std::runtime_error("sfncon: elements.sfncon{4|6|8}s.slave: node index out of range");
               }
          }

          const octave_value k = ov_k.xelem(l);

          if (!(k.is_function_handle() || k.is_inline_function() || k.is_string())) {
               throw std::runtime_error("sfncon: elements.sfncon{4|6|8}s.k must be a function handle, inline function or string");
          }
          // const std::complex<double> sigma0 = ov_sigma0.xelem(l).complex_value();
          // const double sigma_delta = ov_sigma_delta.xelem(l).scalar_value();

          switch (oElemType.type) {
          case ElementTypes::ELEM_SFNCON4S:
               SurfToNodeConstr<ShapeIso4>::BuildContacts(id,
                                                          nodes,
                                                          nidxmaster,
                                                          nidxslave,
                                                          maxdist,
                                                          eConstrType,
                                                          uConstraintFlags,
                                                          eDomain,
                                                          rgElemBlocks,
                                                          k);
               break;
          case ElementTypes::ELEM_SFNCON6S:
               SurfToNodeConstr<ShapeTria6>::BuildContacts(id,
                                                           nodes,
                                                           nidxmaster,
                                                           nidxslave,
                                                           maxdist,
                                                           eConstrType,
                                                           uConstraintFlags,
                                                           eDomain,
                                                           rgElemBlocks,
                                                           k);
               break;
          case ElementTypes::ELEM_SFNCON6HS:
               SurfToNodeConstr<ShapeTria6H>::BuildContacts(id,
                                                            nodes,
                                                            nidxmaster,
                                                            nidxslave,
                                                            maxdist,
                                                            eConstrType,
                                                            uConstraintFlags,
                                                            eDomain,
                                                            rgElemBlocks,
                                                            k);
               break;
          case ElementTypes::ELEM_SFNCON8S:
               SurfToNodeConstr<ShapeQuad8>::BuildContacts(id,
                                                           nodes,
                                                           nidxmaster,
                                                           nidxslave,
                                                           maxdist,
                                                           eConstrType,
                                                           uConstraintFlags,
                                                           eDomain,
                                                           rgElemBlocks,
                                                           k);
               break;
          case ElementTypes::ELEM_SFNCON8RS:
               SurfToNodeConstr<ShapeQuad8r>::BuildContacts(id,
                                                            nodes,
                                                            nidxmaster,
                                                            nidxslave,
                                                            maxdist,
                                                            eConstrType,
                                                            uConstraintFlags,
                                                            eDomain,
                                                            rgElemBlocks,
                                                            k);
               break;
          case ElementTypes::ELEM_SFNCON9S:
               SurfToNodeConstr<ShapeQuad9>::BuildContacts(id,
                                                           nodes,
                                                           nidxmaster,
                                                           nidxslave,
                                                           maxdist,
                                                           eConstrType,
                                                           uConstraintFlags,
                                                           eDomain,
                                                           rgElemBlocks,
                                                           k);
               break;
          case ElementTypes::ELEM_SFNCON10S:
               SurfToNodeConstr<ShapeTria10>::BuildContacts(id,
                                                            nodes,
                                                            nidxmaster,
                                                            nidxslave,
                                                            maxdist,
                                                            eConstrType,
                                                            uConstraintFlags,
                                                            eDomain,
                                                            rgElemBlocks,
                                                            k);
               break;
          default:
               FEM_ASSERT(false);
          }
     }
}

#endif

octave_scalar_map
SurfaceNormalVectorPostProc(const array<bool, ElementTypes::iGetNumTypes()>& rgElemUse,
                            const vector<std::unique_ptr<ElementBlockBase> >& rgElemBlocks,
                            const octave_scalar_map& elements,
                            const Matrix& nodes,
                            PostProcData& oSolution,
                            const ParallelOptions& oParaOpt) {
     constexpr octave_idx_type iNumComp = 3;
     octave_scalar_map mapSurfaceNormalVectorElem;

     const auto iterPartVel = elements.seek("particle_velocity");

     if (iterPartVel == elements.end()) {
          return mapSurfaceNormalVectorElem;
     }

     const octave_value ovPartVel = elements.contents(iterPartVel);

     if (!ovPartVel.isstruct()) {
          throw std::runtime_error("surface normal vector: mesh.elements.particle_velocity must be a scalar struct");
     }

     const octave_scalar_map mapPartVel = ovPartVel.scalar_map_value();

     for (octave_idx_type j = 0; j < ElementTypes::iGetNumTypes(); ++j) {
          const ElementTypes::TypeInfo& oElemType = ElementTypes::GetType(j);

          if (!rgElemUse[oElemType.type]) {
               continue;
          }

          switch (oElemType.type) {
          case ElementTypes::ELEM_PARTICLE_VEL_ISO4:
          case ElementTypes::ELEM_PARTICLE_VEL_QUAD8:
          case ElementTypes::ELEM_PARTICLE_VEL_QUAD8R:
          case ElementTypes::ELEM_PARTICLE_VEL_QUAD9:
          case ElementTypes::ELEM_PARTICLE_VEL_TRIA6:
          case ElementTypes::ELEM_PARTICLE_VEL_TRIA6H:
          case ElementTypes::ELEM_PARTICLE_VEL_TRIA10: {
               const auto iterElem = mapPartVel.seek(oElemType.name);

               if (iterElem == mapPartVel.end()) {
                    continue;
               }

               const octave_value ovElem = mapPartVel.contents(iterElem);

               if (!(ovElem.isstruct() && ovElem.numel() == 1)) {
                    throw std::runtime_error("surface normal vector: mesh.elements.particle_velocity."s + oElemType.name + " must be a scalar struct");
               }

               const octave_scalar_map mapElem = ovElem.scalar_map_value();

               const auto iterNodes = mapElem.seek("nodes");

               if (iterNodes == mapElem.end()) {
                    throw std::runtime_error("surface normal vector: field mesh.elements.particle_velocity."s + oElemType.name + ".nodes not found");
               }
               const octave_value ovNodes = mapElem.contents(iterNodes);

               if (!(ovNodes.is_matrix_type() && ovNodes.isinteger() && ovNodes.columns() == oElemType.min_nodes)) {
                    throw std::runtime_error("surface normal vector: number of columns does not match for field mesh.elements.particle_velocity."s + oElemType.name + ".nodes");
               }

               const int32NDArray elemNodes = ovNodes.int32_array_value();

               oSolution.SetField(PostProcData::VEC_EL_SURFACE_NORMAL_VECTOR_RE,
                                  oElemType.type,
                                  NDArray(dim_vector(elemNodes.rows(),
                                                     elemNodes.columns(),
                                                     iNumComp),
                                          0.));


               for (const auto& pElemBlock: rgElemBlocks) {
                    if (pElemBlock->GetElementType() == oElemType.type) {
                         pElemBlock->PostProcElem(Element::VEC_SURFACE_NORMAL_VECTOR, oSolution, oParaOpt);
                    }
               }

               const NDArray nel = oSolution.GetField(PostProcData::VEC_EL_SURFACE_NORMAL_VECTOR_RE, oElemType.type);

               mapSurfaceNormalVectorElem.assign(oElemType.name, nel);
          } break;

          default:
               break;
          }
     }

     return mapSurfaceNormalVectorElem;
}

octave_scalar_map
SurfaceAreaPostProc(const array<bool, ElementTypes::iGetNumTypes()>& rgElemUse,
                    const vector<std::unique_ptr<ElementBlockBase> >& rgElemBlocks,
                    const Matrix& nodes,
                    PostProcData& oSolution,
                    const ParallelOptions& oParaOpt) {
     octave_scalar_map mapSurfaceAreaElem;

     for (const auto& pElemBlock: rgElemBlocks) {
          const ElementTypes::TypeId eltype = pElemBlock->GetElementType();

          switch (eltype) {
          case ElementTypes::ELEM_PRESSURE_ISO4:
          case ElementTypes::ELEM_PRESSURE_QUAD8:
          case ElementTypes::ELEM_PRESSURE_QUAD8R:
          case ElementTypes::ELEM_PRESSURE_QUAD9:
          case ElementTypes::ELEM_PRESSURE_TRIA6:
          case ElementTypes::ELEM_PRESSURE_TRIA6H:
          case ElementTypes::ELEM_PRESSURE_TRIA10:
               break;

          default:
               continue;
          }

          const ElementTypes::TypeInfo& oElemType = ElementTypes::GetType(eltype);

          FEM_ASSERT(oElemType.type == eltype);
          FEM_ASSERT(oElemType.min_nodes == oElemType.max_nodes);
          FEM_ASSERT(oElemType.min_nodes > 0);

          const octave_idx_type iNumNodes = oElemType.min_nodes;
          const octave_idx_type iNumElem = pElemBlock->iGetNumElem();

          oSolution.SetField(PostProcData::VEC_EL_SURFACE_AREA_RE,
                             oElemType.type,
                             NDArray(dim_vector(iNumElem, iNumNodes), 0.));

          pElemBlock->PostProcElem(Element::VEC_SURFACE_AREA, oSolution, oParaOpt);

          const NDArray Aelem = oSolution.GetField(PostProcData::VEC_EL_SURFACE_AREA_RE, oElemType.type);

          mapSurfaceAreaElem.assign(oElemType.name, Aelem);
     }

     return mapSurfaceAreaElem;
}

template <typename T>
octave_scalar_map AcousticPostProc(const array<bool, ElementTypes::iGetNumTypes()>& rgElemUse,
                                   const vector<std::unique_ptr<ElementBlockBase> >& rgElemBlocks,
                                   const octave_scalar_map& elements,
                                   const Matrix& nodes,
                                   PostProcData& oSolution,
                                   const Element::FemMatrixType eMatType,
                                   const ParallelOptions& oParaOpt) {
     constexpr octave_idx_type iNumComp = 3;
     const octave_idx_type iNumLoads = oSolution.GetNumSteps();
     const octave_idx_type iNumNodes = nodes.rows();
     typedef typename PostProcTypeTraits<T>::NDArrayType TNDArray;
     typedef PostProcFieldHelper<T> PPFH;
     TNDArray oNodalVelocity(dim_vector(iNumNodes, iNumComp, iNumLoads), 0.);
     octave_scalar_map mapVelocityElemVec;

     for (octave_idx_type m = 0; m < ElementTypes::iGetNumTypes(); ++m) {
          const ElementTypes::TypeInfo& oElemType = ElementTypes::GetType(m);

          if (!rgElemUse[oElemType.type]) {
               continue;
          }

          switch (oElemType.type) {
          case ElementTypes::ELEM_ISO8:
          case ElementTypes::ELEM_ISO20:
          case ElementTypes::ELEM_ISO20R:
          case ElementTypes::ELEM_ISO27:
          case ElementTypes::ELEM_PENTA15:
          case ElementTypes::ELEM_TET10H:
          case ElementTypes::ELEM_TET10:
          case ElementTypes::ELEM_TET20: {
               const auto iter_elem = elements.seek(oElemType.name);

               if (iter_elem == elements.end()) {
                    continue;
               }

               const int32NDArray elem_nodes = elements.contents(iter_elem).int32_array_value();

               const octave_idx_type iNumElem = elem_nodes.rows();
               const octave_idx_type iNumNodesElem = elem_nodes.columns();

               oSolution.SetField(PPFH::ConvertFieldType(PostProcData::VEC_EL_ACOUSTIC_PART_VEL_RE),
                                  oElemType.type,
                                  TNDArray(dim_vector(iNumElem,
                                                      iNumNodesElem,
                                                      iNumComp,
                                                      iNumLoads),
                                           T{}));


               for (const auto& pElemBlock: rgElemBlocks) {
                    if (pElemBlock->GetElementType() == oElemType.type) {
                         pElemBlock->PostProcElem(PPFH::ConvertMatrixType(Element::VEC_PARTICLE_VELOCITY), oSolution, oParaOpt);
                    }
               }

               const TNDArray oElemVelocity = oSolution.GetField(PPFH::ConvertFieldType(PostProcData::VEC_EL_ACOUSTIC_PART_VEL_RE), oElemType.type);

               for (octave_idx_type l = 0; l < iNumLoads; ++l) {
                    for (octave_idx_type k = 0; k < iNumComp; ++k) {
                         for (octave_idx_type j = 0; j < iNumNodesElem; ++j) {
                              for (octave_idx_type i = 0; i < iNumElem; ++i) {
                                   const octave_idx_type inode = elem_nodes.xelem(i + iNumElem * j).value() - 1;
                                   oNodalVelocity.xelem(inode + iNumNodes * (k + iNumComp * l)) = oElemVelocity.xelem(i + iNumElem * (j + iNumNodesElem * (l * iNumComp + k)));
                              }
                         }
                    }
               }

               mapVelocityElemVec.assign(oElemType.name, oElemVelocity);
          } break;

          default:
               break;
          }
     }

     octave_scalar_map mapAcoustics;

     mapAcoustics.assign("v", mapVelocityElemVec);

     oSolution.SetField(PPFH::ConvertFieldType(PostProcData::VEC_NO_ACOUSTIC_PART_VEL_RE),
                        ElementTypes::ELEM_TYPE_UNKNOWN,
                        oNodalVelocity);

     const auto iter_bound = elements.seek("acoustic_boundary");

     if (iter_bound != elements.end()) {
          octave_scalar_map mapVelocityElemNorm, mapAcousticIntensityElem, mapSoundPowerElem;
          const octave_scalar_map mapAcousticBoundary = elements.contents(iter_bound).scalar_map_value();

          for (octave_idx_type j = 0; j < ElementTypes::iGetNumTypes(); ++j) {
               const ElementTypes::TypeInfo& oElemType = ElementTypes::GetType(j);

               if (!rgElemUse[oElemType.type]) {
                    continue;
               }

               switch (oElemType.type) {
               case ElementTypes::ELEM_ACOUSTIC_BND_ISO4:
               case ElementTypes::ELEM_ACOUSTIC_BND_QUAD8:
               case ElementTypes::ELEM_ACOUSTIC_BND_QUAD8R:
               case ElementTypes::ELEM_ACOUSTIC_BND_QUAD9:
               case ElementTypes::ELEM_ACOUSTIC_BND_TRIA6:
               case ElementTypes::ELEM_ACOUSTIC_BND_TRIA6H:
               case ElementTypes::ELEM_ACOUSTIC_BND_TRIA10: {
                    const auto iter_elem = mapAcousticBoundary.seek(oElemType.name);

                    if (iter_elem == mapAcousticBoundary.end()) {
                         continue;
                    }

                    const int32NDArray elem_nodes = mapAcousticBoundary.contents(iter_elem).int32_array_value();

                    const octave_idx_type iNumElem = elem_nodes.rows();
                    const octave_idx_type iNumNodesElem = elem_nodes.columns();

                    switch (eMatType) {
                    case Element::SCA_ACOUSTIC_INTENSITY:
                    case Element::SCA_ACOUSTIC_INTENSITY_C:
                         oSolution.SetField(PostProcData::SCA_EL_ACOUSTIC_INTENSITY_RE,
                                            oElemType.type,
                                            NDArray(dim_vector(iNumElem,
                                                               iNumNodesElem,
                                                               iNumLoads),
                                                    0.));
                         oSolution.SetField(PostProcData::SCA_EL_ACOUSTIC_SOUND_POWER_RE,
                                            oElemType.type,
                                            NDArray(dim_vector(iNumElem,
                                                               iNumLoads),
                                                    0.));
                         break;
                    case Element::VEC_PARTICLE_VELOCITY:
                    case Element::VEC_PARTICLE_VELOCITY_C:
                         oSolution.SetField(PPFH::ConvertFieldType(PostProcData::SCA_EL_ACOUSTIC_PART_VEL_NORM_RE),
                                            oElemType.type,
                                            TNDArray(dim_vector(iNumElem,
                                                                iNumNodesElem,
                                                                iNumLoads),
                                                     T{}));
                         break;
                    default:
                         FEM_ASSERT(0);
                    }

                    for (const auto& pElemBlock: rgElemBlocks) {
                         if (pElemBlock->GetElementType() == oElemType.type) {
                              pElemBlock->PostProcElem(eMatType, oSolution, oParaOpt);
                         }
                    }

                    switch (eMatType) {
                    case Element::SCA_ACOUSTIC_INTENSITY:
                    case Element::SCA_ACOUSTIC_INTENSITY_C:
                         mapAcousticIntensityElem.assign(oElemType.name,
                                                         oSolution.GetField(PostProcData::SCA_EL_ACOUSTIC_INTENSITY_RE,
                                                                            oElemType.type));
                         mapSoundPowerElem.assign(oElemType.name,
                                                  oSolution.GetField(PostProcData::SCA_EL_ACOUSTIC_SOUND_POWER_RE,
                                                                     oElemType.type));
                         break;
                    case Element::VEC_PARTICLE_VELOCITY:
                    case Element::VEC_PARTICLE_VELOCITY_C:
                         mapVelocityElemNorm.assign(oElemType.name,
                                                    oSolution.GetField(PPFH::ConvertFieldType(PostProcData::SCA_EL_ACOUSTIC_PART_VEL_NORM_RE),
                                                                       oElemType.type));
                         break;
                    default:
                         FEM_ASSERT(0);
                    };

               } break;

               default:
                    break;
               }
          }

          switch (eMatType) {
          case Element::SCA_ACOUSTIC_INTENSITY:
          case Element::SCA_ACOUSTIC_INTENSITY_C:
               mapAcoustics.assign("I", mapAcousticIntensityElem);
               mapAcoustics.assign("P", mapSoundPowerElem);
               break;
          case Element::VEC_PARTICLE_VELOCITY:
          case Element::VEC_PARTICLE_VELOCITY_C:
               mapAcoustics.assign("vn", mapVelocityElemNorm);
               break;
          default:
               FEM_ASSERT(0);
          }
     }

     return mapAcoustics;
}

// PKG_ADD: autoload("fem_ass_matrix", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("fem_ass_dof_map", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("fem_pre_mesh_constr_surf_to_node", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_MAT_MASS", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_MAT_MASS_SYM", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_MAT_MASS_SYM_L", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_MAT_MASS_LUMPED", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_MAT_STIFFNESS", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_MAT_STIFFNESS_SYM", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_MAT_STIFFNESS_SYM_L", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_MAT_STIFFNESS_IM", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_MAT_STIFFNESS_TAU0", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_MAT_STIFFNESS_OMEGA", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_MAT_STIFFNESS_OMEGA_DOT", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_MAT_DAMPING_OMEGA", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_MAT_DAMPING", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_MAT_DAMPING_SYM", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_MAT_DAMPING_SYM_L", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_SCA_TOT_MASS", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_VEC_INERTIA_M1", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_MAT_INERTIA_J", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_MAT_INERTIA_INV3", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_MAT_INERTIA_INV4", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_MAT_INERTIA_INV5", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_MAT_INERTIA_INV8", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_MAT_INERTIA_INV9", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_VEC_LOAD_CONSISTENT", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_VEC_LOAD_LUMPED", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_VEC_STRESS_CAUCH", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_VEC_STRAIN_TOTAL", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_SCA_STRESS_VMIS", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_CT_FIXED", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_CT_SLIDING", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_MAT_THERMAL_COND", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_MAT_HEAT_CAPACITY", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_VEC_LOAD_THERMAL", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_MAT_MASS_ACOUSTICS_RE", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_MAT_MASS_ACOUSTICS_IM", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_MAT_STIFFNESS_ACOUSTICS_RE", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_MAT_STIFFNESS_ACOUSTICS_IM", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_MAT_DAMPING_ACOUSTICS_RE", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_MAT_DAMPING_ACOUSTICS_IM", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_VEC_LOAD_ACOUSTICS", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_VEC_PARTICLE_VELOCITY", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_VEC_PARTICLE_VELOCITY_C", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_SCA_ACOUSTIC_INTENSITY", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_SCA_ACOUSTIC_INTENSITY_C", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_VEC_SURFACE_NORMAL_VECTOR", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_VEC_SURFACE_AREA", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_MAT_MASS_FLUID_STRUCT_RE", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_MAT_MASS_FLUID_STRUCT_IM", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_MAT_STIFFNESS_FLUID_STRUCT_RE", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_MAT_STIFFNESS_FLUID_STRUCT_IM", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_MAT_DAMPING_FLUID_STRUCT_RE", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_MAT_DAMPING_FLUID_STRUCT_IM", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_VEC_LOAD_FLUID_STRUCT", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_DO_THERMAL", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_DO_STRUCTURAL", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_DO_ACOUSTICS", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_DO_FLUID_STRUCT", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_VEC_COLL_MASS", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_VEC_COLL_STIFFNESS", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_VEC_COLL_HEAT_CAPACITY", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_VEC_COLL_THERMAL_COND", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_VEC_COLL_MASS_ACOUSTICS", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_VEC_COLL_STIFF_ACOUSTICS", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_VEC_COLL_MASS_FLUID_STRUCT", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_VEC_COLL_STIFF_FLUID_STRUCT", "__mboct_fem_pkg__.oct");

// PKG_DEL: autoload("fem_ass_matrix", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("fem_ass_dof_map", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("fem_pre_mesh_constr_surf_to_node", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_MAT_MASS", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_MAT_MASS_SYM", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_MAT_MASS_SYM_L", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_MAT_MASS_LUMPED", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_MAT_STIFFNESS", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_MAT_STIFFNESS_SYM", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_MAT_STIFFNESS_SYM_L", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_MAT_STIFFNESS_TAU0", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_MAT_STIFFNESS_OMEGA", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_MAT_STIFFNESS_OMEGA_DOT", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_MAT_DAMPING_OMEGA", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_MAT_DAMPING", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_MAT_DAMPING_SYM", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_MAT_DAMPING_SYM_L", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_SCA_TOT_MASS", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_VEC_INERTIA_M1", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_MAT_INERTIA_J", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_MAT_INERTIA_INV3", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_MAT_INERTIA_INV4", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_MAT_INERTIA_INV5", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_MAT_INERTIA_INV8", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_MAT_INERTIA_INV9", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_VEC_LOAD_CONSISTENT", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_VEC_LOAD_LUMPED", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_VEC_STRESS_CAUCH", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_VEC_STRAIN_TOTAL", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_SCA_STRESS_VMIS", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_CT_FIXED", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_CT_SLIDING", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_MAT_THERMAL_COND", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_MAT_HEAT_CAPACITY", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_VEC_LOAD_THERMAL", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_MAT_MASS_ACOUSTICS_RE", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_MAT_MASS_ACOUSTICS_IM", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_MAT_STIFFNESS_ACOUSTICS_RE", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_MAT_STIFFNESS_ACOUSTICS_IM", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_MAT_DAMPING_ACOUSTICS_RE", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_MAT_DAMPING_ACOUSTICS_IM", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_VEC_LOAD_ACOUSTICS", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_VEC_PARTICLE_VELOCITY", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_VEC_PARTICLE_VELOCITY_C", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_SCA_ACOUSTIC_INTENSITY", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_SCA_ACOUSTIC_INTENSITY_C", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_VEC_SURFACE_NORMAL_VECTOR", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_VEC_SURFACE_AREA", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_MAT_MASS_FLUID_STRUCT_RE", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_MAT_MASS_FLUID_STRUCT_IM", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_MAT_STIFFNESS_FLUID_STRUCT_RE", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_MAT_STIFFNESS_FLUID_STRUCT_IM", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_MAT_DAMPING_FLUID_STRUCT_RE", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_MAT_DAMPING_FLUID_STRUCT_IM", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_VEC_LOAD_FLUID_STRUCT", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_DO_THERMAL", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_DO_STRUCTURAL", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_DO_ACOUSTICS", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_DO_FLUID_STRUCT", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_VEC_COLL_MASS", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_VEC_COLL_STIFFNESS", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_VEC_COLL_HEAT_CAPACITY", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_VEC_COLL_THERMAL_COND", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_VEC_COLL_MASS_ACOUSTICS", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_VEC_COLL_STIFF_ACOUSTICS", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_VEC_COLL_MASS_FLUID_STRUCT", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_VEC_COLL_STIFF_FLUID_STRUCT", "__mboct_fem_pkg__.oct", "remove");

DEFUN_DLD(fem_ass_dof_map, args, nargout,
          "-*- texinfo -*-\n"
          "@deftypefn {} @var{dof_map} = fem_ass_dof_map(@var{mesh}, @var{load_case})\n"
          "Build a data structure which assigns global degrees of freedom to each node and each element with Lagrange multipliers.\n\n"
          "@var{mesh} @dots{} Finite element mesh data structure\n\n"
          "@var{load_case} @dots{} Scalar struct for defining nodal constraints.\n\n"
          "@var{dof_map} @dots{} Degree of freedom mapping\n\n"
          "@seealso{fem_tests}\n\n"
          "@end deftypefn\n")
{
     octave_value_list retval;

     try {
          if (args.length() < 2) {
               print_usage();
               return retval;
          }

          const octave_scalar_map m_mesh = args(0).scalar_map_value();

#if OCTAVE_MAJOR_VERSION < 6
          if (error_state) {
               return retval;
          }
#endif

          const auto it_nodes = m_mesh.seek("nodes");

          if (it_nodes == m_mesh.end()) {
               error("field mesh.nodes not found");
               return retval;
          }

          const Matrix nodes = m_mesh.contents(it_nodes).matrix_value();

#if OCTAVE_MAJOR_VERSION < 6
          if (error_state) {
               return retval;
          }
#endif

          if (nodes.ndims() != 2 || nodes.columns() != 6) {
               error("mesh.nodes must be a Nx6 matrix");
               return retval;
          }

          const octave_idx_type iNumNodes = nodes.rows();

          const octave_scalar_map m_load_case = args(1).scalar_map_value();

#if OCTAVE_MAJOR_VERSION < 6
          if (error_state) {
               return retval;
          }
#endif

          const auto it_locked_dof = m_load_case.seek("locked_dof");

          if (it_locked_dof == m_load_case.end()) {
               error("field load_case.locked_dof not found");
               return retval;
          }

          const boolNDArray locked_dof = m_load_case.contents(it_locked_dof).bool_array_value();

#if OCTAVE_MAJOR_VERSION < 6
          if (error_state) {
               return retval;
          }
#endif
          DofMap::DomainType eDomain = DofMap::DO_STRUCTURAL;

          const auto it_domain = m_load_case.seek("domain");

          if (it_domain != m_load_case.end()) {
               eDomain = static_cast<DofMap::DomainType>(m_load_case.contents(it_domain).int_value());
          }

          const octave_idx_type iNodeMaxDofIndex = DofMap::iGetNodeMaxDofIndex(eDomain);

          if (iNodeMaxDofIndex <= 0) {
               error("invalid value for load_case.domain=%d", static_cast<int>(eDomain));
               return retval;
          }

          if (locked_dof.ndims() != 2 ||
              locked_dof.columns() != iNodeMaxDofIndex ||
              locked_dof.rows() != iNumNodes) {
               error("size of load_case.locked_dof is not valid");
               return retval;
          }

          const auto iter_elements = m_mesh.seek("elements");

          if (iter_elements == m_mesh.end()) {
               error("field mesh.elements not found");
               return retval;
          }

          const octave_scalar_map m_elements = m_mesh.contents(iter_elements).scalar_map_value();

#if OCTAVE_MAJOR_VERSION < 6
          if (error_state) {
               return retval;
          }
#endif

          const auto it_materials = m_mesh.seek("materials");

          if (it_materials == m_mesh.end()) {
               throw std::runtime_error("fem_ass_dof_map: missing field mesh.materials in argument mesh");
          }

          const octave_scalar_map materials(m_mesh.contents(it_materials).scalar_map_value());

#if OCTAVE_MAJOR_VERSION < 6
          if (error_state) {
               throw std::runtime_error("fem_ass_dof_map: mesh.materials must be a scalar struct in argument mesh");
          }
#endif

          const auto it_material_data = m_mesh.seek("material_data");

          if (it_material_data == m_mesh.end()) {
               throw std::runtime_error("fem_ass_dof_map: missing field mesh.material_data in argument mesh");
          }

          const octave_map material_data(m_mesh.contents(it_material_data).map_value());

#if OCTAVE_MAJOR_VERSION < 6
          if (error_state) {
               throw std::runtime_error("fem_ass_dof_map: mesh.material_data must be a struct array in argument mesh");
          }
#endif

          const vector<Material> rgMaterials = Material::ExtractMaterialData(material_data, eDomain);

          const auto it_dof_in_use = m_load_case.seek("dof_in_use");

          boolNDArray dof_in_use(it_dof_in_use == m_load_case.end()
                                 ? boolNDArray(dim_vector(iNumNodes, iNodeMaxDofIndex), false)
                                 : m_load_case.contents(it_dof_in_use).bool_array_value());

          if (dof_in_use.ndims() != 2 || dof_in_use.rows() != iNumNodes || dof_in_use.columns() != iNodeMaxDofIndex) {
               throw std::runtime_error("fem_ass_dof_map: invalid size for dof_in_use");
          }

          for (octave_idx_type i = 0; i < ElementTypes::iGetNumTypes(); ++i) {
               const auto& oElemType = ElementTypes::GetType(i);

               const auto iter_elem_type = m_elements.seek(oElemType.name);

               if (iter_elem_type == m_elements.end()) {
                    continue;
               }

               switch (oElemType.type) {
               case ElementTypes::ELEM_ISO8:
               case ElementTypes::ELEM_ISO20:
               case ElementTypes::ELEM_ISO20R:
               case ElementTypes::ELEM_ISO27:
               case ElementTypes::ELEM_PENTA15:
               case ElementTypes::ELEM_TET10H:
               case ElementTypes::ELEM_TET10:
               case ElementTypes::ELEM_TET20: {
                    const auto iter_elem_mat = materials.seek(oElemType.name);

                    if (iter_elem_mat == materials.end()) {
                         throw std::runtime_error("fem_ass_dof_map: missing field mesh.materials."s + oElemType.name + " in argument mesh");
                    }

                    const int32NDArray elem_mat = materials.contents(iter_elem_mat).int32_array_value();

                    if (elem_mat.columns() != 1) {
                         throw std::runtime_error("fem_ass_dof_map: invalid number of columns in mesh.materials."s + oElemType.name);
                    }

                    const int32NDArray elnodes = m_elements.contents(iter_elem_type).int32_array_value();

                    if (!(elnodes.columns() >= oElemType.min_nodes && elnodes.columns() <= oElemType.max_nodes)) {
                         throw std::runtime_error("fem_ass_dof_map: invalid number of columns in mesh.elements."s + oElemType.name);
                    }

                    const octave_idx_type iNumElem = elnodes.rows();
                    const octave_idx_type iNumNodesElem = elnodes.columns();

                    if (elem_mat.rows() != iNumElem) {
                         throw std::runtime_error("fem_ass_dof_map: inconsistent size of mesh.elements."s + oElemType.name + " and mesh.materials." + oElemType.name);
                    }

                    for (octave_idx_type j = 0; j < iNumElem; ++j) {
                         const size_t imaterial = elem_mat.xelem(j).value() - 1;

                         if (imaterial >= rgMaterials.size()) {
                              throw std::runtime_error("fem_ass_dof_map: invalid index in field mesh.materials."s + oElemType.name);
                         }

                         const Material::MatType eMatType = rgMaterials[imaterial].GetMaterialType();

                         for (octave_idx_type k = 0; k < iNumNodesElem; ++k) {
                              const octave_idx_type idxnode = elnodes.xelem(j + iNumElem * k).value();

                              if (idxnode < 1 || idxnode > iNumNodes) {
                                   error("node index %Ld of element mesh.elements.%s(%Ld, %Ld) out of range %Ld:%Ld",
                                         static_cast<long long>(idxnode),
                                         oElemType.name,
                                         static_cast<long long>(j),
                                         static_cast<long long>(k),
                                         1LL,
                                         static_cast<long long>(iNumNodes));
                                   return retval;
                              }

                              octave_idx_type iNodeDofMin = 1, iNodeDofMax = -1;

                              switch (eDomain) {
                              case DofMap::DO_STRUCTURAL:
                                   iNodeDofMin = 0;
                                   iNodeDofMax = 2;
                                   break;
                              case DofMap::DO_THERMAL:
                              case DofMap::DO_ACOUSTICS:
                                   iNodeDofMin = iNodeDofMax = 0;
                                   break;
                              case DofMap::DO_FLUID_STRUCT:
                                   switch (eMatType) {
                                   case Material::MAT_TYPE_SOLID:
                                        iNodeDofMin = 0;
                                        iNodeDofMax = 2;
                                        break;
                                   case Material::MAT_TYPE_FLUID:
                                        iNodeDofMin = iNodeDofMax = 6;
                                        break;
                                   default:
                                        throw std::logic_error("fem_ass_dof_map: invalid material type for fluid structure interaction");
                                   }
                                   break;
                              default:
                                   throw std::runtime_error("fem_ass_dof_map: unknown value for dof_map.domain");
                              }

                              for (octave_idx_type l = iNodeDofMin; l <= iNodeDofMax; ++l) {
                                   dof_in_use.xelem(idxnode - 1, l) = true;
                              }
                         }
                    }
               } break;
               case ElementTypes::ELEM_THERM_CONV_ISO4:
               case ElementTypes::ELEM_THERM_CONV_QUAD8:
               case ElementTypes::ELEM_THERM_CONV_QUAD8R:
               case ElementTypes::ELEM_THERM_CONV_QUAD9:
               case ElementTypes::ELEM_THERM_CONV_TRIA6:
               case ElementTypes::ELEM_THERM_CONV_TRIA6H:
               case ElementTypes::ELEM_THERM_CONV_TRIA10:
               case ElementTypes::ELEM_PARTICLE_VEL_ISO4:
               case ElementTypes::ELEM_PARTICLE_VEL_QUAD8:
               case ElementTypes::ELEM_PARTICLE_VEL_QUAD8R:
               case ElementTypes::ELEM_PARTICLE_VEL_QUAD9:
               case ElementTypes::ELEM_PARTICLE_VEL_TRIA6:
               case ElementTypes::ELEM_PARTICLE_VEL_TRIA6H:
               case ElementTypes::ELEM_PARTICLE_VEL_TRIA10:
               case ElementTypes::ELEM_ACOUSTIC_IMPE_ISO4:
               case ElementTypes::ELEM_ACOUSTIC_IMPE_QUAD8:
               case ElementTypes::ELEM_ACOUSTIC_IMPE_QUAD8R:
               case ElementTypes::ELEM_ACOUSTIC_IMPE_QUAD9:
               case ElementTypes::ELEM_ACOUSTIC_IMPE_TRIA6:
               case ElementTypes::ELEM_ACOUSTIC_IMPE_TRIA6H:
               case ElementTypes::ELEM_ACOUSTIC_IMPE_TRIA10:
               case ElementTypes::ELEM_FLUID_STRUCT_ISO4:
               case ElementTypes::ELEM_FLUID_STRUCT_QUAD8:
               case ElementTypes::ELEM_FLUID_STRUCT_QUAD8R:
               case ElementTypes::ELEM_FLUID_STRUCT_QUAD9:
               case ElementTypes::ELEM_FLUID_STRUCT_TRIA6:
               case ElementTypes::ELEM_FLUID_STRUCT_TRIA6H:
               case ElementTypes::ELEM_FLUID_STRUCT_TRIA10:
                    if (eDomain == DofMap::DO_THERMAL || eDomain == DofMap::DO_ACOUSTICS || eDomain == DofMap::DO_FLUID_STRUCT) {
                         static constexpr char elem_name[][23] = {"convection", "particle_velocity", "acoustic_impedance", "fluid_struct_interface"};

                         enum {
                              EL_IDX_CONVECTION,
                              EL_IDX_PARTICLE_VEL,
                              EL_IDX_ACOUSTIC_IMPEDANCE,
                              EL_IDX_FLUID_STRUCT_INTERFACE
                         } ielem_name;

                         switch (oElemType.type) {
                         case ElementTypes::ELEM_THERM_CONV_ISO4:
                         case ElementTypes::ELEM_THERM_CONV_QUAD8:
                         case ElementTypes::ELEM_THERM_CONV_QUAD8R:
                         case ElementTypes::ELEM_THERM_CONV_QUAD9:
                         case ElementTypes::ELEM_THERM_CONV_TRIA6:
                         case ElementTypes::ELEM_THERM_CONV_TRIA6H:
                         case ElementTypes::ELEM_THERM_CONV_TRIA10:
                              ielem_name = EL_IDX_CONVECTION;
                              break;

                         case ElementTypes::ELEM_PARTICLE_VEL_ISO4:
                         case ElementTypes::ELEM_PARTICLE_VEL_QUAD8:
                         case ElementTypes::ELEM_PARTICLE_VEL_QUAD8R:
                         case ElementTypes::ELEM_PARTICLE_VEL_QUAD9:
                         case ElementTypes::ELEM_PARTICLE_VEL_TRIA6:
                         case ElementTypes::ELEM_PARTICLE_VEL_TRIA6H:
                         case ElementTypes::ELEM_PARTICLE_VEL_TRIA10:
                              ielem_name = EL_IDX_PARTICLE_VEL;
                              break;

                         case ElementTypes::ELEM_ACOUSTIC_IMPE_ISO4:
                         case ElementTypes::ELEM_ACOUSTIC_IMPE_QUAD8:
                         case ElementTypes::ELEM_ACOUSTIC_IMPE_QUAD8R:
                         case ElementTypes::ELEM_ACOUSTIC_IMPE_QUAD9:
                         case ElementTypes::ELEM_ACOUSTIC_IMPE_TRIA6:
                         case ElementTypes::ELEM_ACOUSTIC_IMPE_TRIA6H:
                         case ElementTypes::ELEM_ACOUSTIC_IMPE_TRIA10:
                              ielem_name = EL_IDX_ACOUSTIC_IMPEDANCE;
                              break;
                         case ElementTypes::ELEM_FLUID_STRUCT_ISO4:
                         case ElementTypes::ELEM_FLUID_STRUCT_QUAD8:
                         case ElementTypes::ELEM_FLUID_STRUCT_QUAD8R:
                         case ElementTypes::ELEM_FLUID_STRUCT_QUAD9:
                         case ElementTypes::ELEM_FLUID_STRUCT_TRIA6:
                         case ElementTypes::ELEM_FLUID_STRUCT_TRIA6H:
                         case ElementTypes::ELEM_FLUID_STRUCT_TRIA10:
                              ielem_name = EL_IDX_FLUID_STRUCT_INTERFACE;
                              break;
                         default:
                              FEM_ASSERT(0);
                              throw std::logic_error("fem_ass_dof_map: unexpected element type");
                         }

                         const auto iter_elem_name = m_elements.seek(elem_name[ielem_name]);

                         if (iter_elem_name == m_elements.end()) {
                              break;
                         }

                         const octave_value ov_elem_name = m_elements.contents(iter_elem_name);

                         if (!(ov_elem_name.isstruct() && ov_elem_name.numel() == 1)) {
                              throw std::runtime_error("fem_ass_dof_map: mesh.elements."s + elem_name[ielem_name] + " must be a scalar struct");
                         }

                         const octave_scalar_map m_elem_name = ov_elem_name.scalar_map_value();

                         const auto iter_surf_elem_type = m_elem_name.seek(oElemType.name);

                         if (iter_surf_elem_type == m_elem_name.end()) {
                              break;
                         }

                         const octave_value ov_elem_type = m_elem_name.contents(iter_surf_elem_type);

                         octave_value ov_elnodes;

                         switch (oElemType.type) {
                         case ElementTypes::ELEM_FLUID_STRUCT_ISO4:
                         case ElementTypes::ELEM_FLUID_STRUCT_QUAD8:
                         case ElementTypes::ELEM_FLUID_STRUCT_QUAD8R:
                         case ElementTypes::ELEM_FLUID_STRUCT_QUAD9:
                         case ElementTypes::ELEM_FLUID_STRUCT_TRIA6:
                         case ElementTypes::ELEM_FLUID_STRUCT_TRIA6H:
                         case ElementTypes::ELEM_FLUID_STRUCT_TRIA10: {
                              ov_elnodes = ov_elem_type;

                              if (!(ov_elnodes.is_matrix_type() && ov_elnodes.isinteger() && ov_elnodes.columns() == oElemType.max_nodes)) {
                                   throw std::runtime_error("fem_ass_dof_map: mesh.elements."s + elem_name[ielem_name] + "." + oElemType.name + " must be an integer matrix");
                              }
                         } break;
                         default: {
                              if (!(ov_elem_type.isstruct() && ov_elem_type.numel() == 1)) {
                                   throw std::runtime_error("fem_ass_dof_map: mesh.elements."s + elem_name[ielem_name] + "." + oElemType.name + " must be a scalar struct");
                              }

                              const octave_scalar_map m_elem_type = ov_elem_type.scalar_map_value();

                              const auto iter_elnodes = m_elem_type.seek("nodes");

                              if (iter_elnodes == m_elem_type.end()) {
                                   throw std::runtime_error("fem_ass_dof_map: missing field mesh.elements."s + elem_name[ielem_name] + "." + oElemType.name + ".nodes");
                              }

                              ov_elnodes = m_elem_type.contents(iter_elnodes);

                              if (!(ov_elnodes.is_matrix_type() && ov_elnodes.isinteger() && ov_elnodes.columns() == oElemType.max_nodes)) {
                                   throw std::runtime_error("fem_ass_dof_map: mesh.elements."s + elem_name[ielem_name] + "." + oElemType.name + ".nodes must be an integer matrix");
                              }
                         }
                         }

                         const int32NDArray elnodes = ov_elnodes.int32_array_value();

                         const octave_idx_type iNumElem = elnodes.rows();
                         const octave_idx_type iNumNodesElem = elnodes.columns();

                         for (octave_idx_type k = 0; k < iNumNodesElem; ++k) {
                              for (octave_idx_type j = 0; j < iNumElem; ++j) {
                                   const octave_idx_type idxnode = elnodes.xelem(j + iNumElem * k).value();

                                   if (idxnode < 1 || idxnode > iNumNodes) {
                                        error("node index %Ld of element mesh.elements.%s.%s(%Ld, %Ld) out of range %Ld:%Ld",
                                              static_cast<long long>(idxnode),
                                              elem_name[ielem_name],
                                              oElemType.name,
                                              static_cast<long long>(j),
                                              static_cast<long long>(k),
                                              1LL,
                                              static_cast<long long>(iNumNodes));
                                        return retval;
                                   }


                                   octave_idx_type iDofIndex = -1;

                                   switch (eDomain) {
                                   case DofMap::DO_ACOUSTICS:
                                   case DofMap::DO_THERMAL:
                                        iDofIndex = 0;
                                        break;
                                   case DofMap::DO_FLUID_STRUCT:
                                        switch (oElemType.type) {
                                        case ElementTypes::ELEM_FLUID_STRUCT_ISO4:
                                        case ElementTypes::ELEM_FLUID_STRUCT_QUAD8:
                                        case ElementTypes::ELEM_FLUID_STRUCT_QUAD8R:
                                        case ElementTypes::ELEM_FLUID_STRUCT_QUAD9:
                                        case ElementTypes::ELEM_FLUID_STRUCT_TRIA6:
                                        case ElementTypes::ELEM_FLUID_STRUCT_TRIA6H:
                                        case ElementTypes::ELEM_FLUID_STRUCT_TRIA10:
                                             for (octave_idx_type l = 0; l < 3; ++l) {
                                                  dof_in_use.xelem(idxnode - 1, l) = true;
                                             }
                                             [[fallthrough]];
                                        default:
                                             iDofIndex = 6;
                                        }
                                        break;
                                   default:
                                        throw std::logic_error("fem_ass_dof_map: domain not supported");
                                   }

                                   dof_in_use.xelem(idxnode - 1, iDofIndex) = true;
                              }
                         }
                    } break;
               case ElementTypes::ELEM_BEAM2:
               case ElementTypes::ELEM_RIGID_BODY:
                    if (eDomain == DofMap::DO_STRUCTURAL || eDomain == DofMap::DO_FLUID_STRUCT) {
                         const octave_map m_elem = m_elements.contents(iter_elem_type).map_value();

#if OCTAVE_MAJOR_VERSION < 6
                         if (error_state) {
                              return retval;
                         }
#endif

                         const auto iter_nodes = m_elem.seek("nodes");

                         if (iter_nodes == m_elem.end()) {
                              error("missing field mesh.elements.%s.nodes", oElemType.name);
                              return retval;
                         }

                         const Cell& ov_nodes = m_elem.contents(iter_nodes);

#if OCTAVE_MAJOR_VERSION < 6
                         if (error_state) {
                              return retval;
                         }
#endif

                         for (octave_idx_type j = 0; j < ov_nodes.numel(); ++j) {
                              const int32NDArray elnodes = ov_nodes(j).int32_array_value();

#if OCTAVE_MAJOR_VERSION < 6
                              if (error_state) {
                                   return retval;
                              }
#endif

                              for (octave_idx_type l = 0; l < elnodes.numel(); ++l) {
                                   const octave_idx_type idxnode = elnodes(l).value();

                                   if (idxnode < 1 || idxnode > iNumNodes) {
                                        error("invalid node index for mesh.elements.%s(%Ld).nodes(%Ld)=%Ld",
                                              oElemType.name,
                                              static_cast<long long>(j),
                                              static_cast<long long>(l),
                                              static_cast<long long>(idxnode));
                                        return retval;
                                   }

                                   for (octave_idx_type k = 0; k < 6; ++k) {
                                        dof_in_use.xelem(idxnode - 1, k) = true;
                                   }
                              }
                         }
                    } break;
               case ElementTypes::ELEM_RBE3:
                    if (eDomain == DofMap::DO_STRUCTURAL || eDomain == DofMap::DO_FLUID_STRUCT) {
                         const octave_map m_rbe3 = m_elements.contents(iter_elem_type).map_value();

#if OCTAVE_MAJOR_VERSION < 6
                         if (error_state) {
                              return retval;
                         }
#endif

                         const auto iter_nodesr = m_rbe3.seek("nodes");

                         if (iter_nodesr == m_rbe3.end()) {
                              error("missing field mesh.elements.rbe3.nodes");
                              return retval;
                         }

                         const Cell& ov_nodesr = m_rbe3.contents(iter_nodesr);

                         for (octave_idx_type j = 0; j < ov_nodesr.numel(); ++j) {
                              const int32NDArray elnodes = ov_nodesr(j).int32_array_value();

#if OCTAVE_MAJOR_VERSION < 6
                              if (error_state) {
                                   return retval;
                              }
#endif

                              if (!elnodes.numel()) {
                                   error("invalid number of nodes for mesh.elements.%s(%Ld).nodes", oElemType.name, static_cast<long long>(j));
                                   return retval;
                              }

                              const octave_idx_type idxnode = elnodes(0).value();

                              if (idxnode < 1 || idxnode > iNumNodes) {
                                   error("invalid node index for mesh.elements.%s(%Ld).nodes(1)=%Ld",
                                         oElemType.name,
                                         static_cast<long long>(j),
                                         static_cast<long long>(idxnode));
                                   return retval;
                              }

                              for (octave_idx_type k = 0; k < 6; ++k) {
                                   dof_in_use.xelem(idxnode - 1, k) = true;
                              }
                         }
                    } break;
               default:
                    continue;
               }
          }

          octave_idx_type icurrdof = 0;

          int32NDArray ndof(dim_vector(iNumNodes, iNodeMaxDofIndex), -1);

          for (octave_idx_type i = 0; i < iNumNodes; ++i) {
               for (octave_idx_type j = 0; j < iNodeMaxDofIndex; ++j) {
                    if (dof_in_use.xelem(i + iNumNodes * j)) {
                         ndof.xelem(i + iNumNodes * j) = locked_dof.xelem(i + iNumNodes * j) ? 0 : ++icurrdof;
                    }
               }
          }

          octave_scalar_map dof_map;

          dof_map.assign("ndof", ndof);
          dof_map.assign("domain", octave_int32{eDomain});

          int32NDArray idx_node(dim_vector(icurrdof, 1), -1);
          octave_idx_type icurrndof = 0;

          for (octave_idx_type i = 0; i < iNumNodes; ++i) {
               for (octave_idx_type j = 0; j < iNodeMaxDofIndex; ++j) {
                    if (ndof.xelem(i + iNumNodes * j).value() > 0) {
                         idx_node.xelem(icurrndof++) = ndof.xelem(i + iNumNodes * j);
                    }
               }
          }

          FEM_ASSERT(icurrndof == icurrdof);

          dof_map.assign("idx_node", idx_node);

          octave_scalar_map m_edof;

          enum {
               CS_JOINT = 0,
               CS_RBE3,
               CS_COUNT
          };

          struct ElemDof {
               int32NDArray dof;
               octave_idx_type elem_count = 0;
               octave_idx_type elem_idx = 0;
               octave_idx_type maxdof = 0;
          };

          array<ElemDof, CS_COUNT> edof;
          octave_idx_type inumlambda = 0;

          enum {
               STAGE_RESERVE = 0,
               STAGE_FILL
          };

          for (octave_idx_type s = STAGE_RESERVE; s <= STAGE_FILL; ++s) {
               for (octave_idx_type i = 0; i < ElementTypes::iGetNumTypes(); ++i) {
                    const auto& oElemType = ElementTypes::GetType(i);

                    switch (oElemType.type) {
                    case ElementTypes::ELEM_RBE3:
                    case ElementTypes::ELEM_JOINT:
                    case ElementTypes::ELEM_SFNCON4:
                    case ElementTypes::ELEM_SFNCON6:
                    case ElementTypes::ELEM_SFNCON6H:
                    case ElementTypes::ELEM_SFNCON8:
                    case ElementTypes::ELEM_SFNCON8R:
                    case ElementTypes::ELEM_SFNCON10:
                    case ElementTypes::ELEM_THERM_CONSTR:
                    case ElementTypes::ELEM_ACOUSTIC_CONSTR:
                         break;

                    default:
                         continue;
                    }

                    const auto iter_etype = m_elements.seek(oElemType.name);

                    if (iter_etype == m_elements.end()) {
                         continue;
                    }

                    const octave_map m_etype = m_elements.contents(iter_etype).map_value();

#if OCTAVE_MAJOR_VERSION < 6
                    if (error_state) {
                         return retval;
                    }
#endif

                    if (!m_etype.numel()) {
                         continue;
                    }

                    switch (oElemType.type) {
                    case ElementTypes::ELEM_JOINT:
                    case ElementTypes::ELEM_SFNCON4:
                    case ElementTypes::ELEM_SFNCON6:
                    case ElementTypes::ELEM_SFNCON6H:
                    case ElementTypes::ELEM_SFNCON8:
                    case ElementTypes::ELEM_SFNCON8R:
                    case ElementTypes::ELEM_SFNCON10:
                    case ElementTypes::ELEM_THERM_CONSTR:
                    case ElementTypes::ELEM_ACOUSTIC_CONSTR: {
                         const char* const fn = (oElemType.type == ElementTypes::ELEM_JOINT ||
                                                 oElemType.type == ElementTypes::ELEM_THERM_CONSTR ||
                                                 oElemType.type == ElementTypes::ELEM_ACOUSTIC_CONSTR)
                              ? "C"
                              : "slave";

                         const auto iter_fn = m_etype.seek(fn);

                         if (iter_fn == m_etype.end()) {
                              error("missing field mesh.elements.%s.%s", oElemType.name, fn);
                              return retval;
                         }

                         const Cell ov_fn = m_etype.contents(iter_fn);

                         Cell ov_constr;

                         switch (oElemType.type) {
                         case ElementTypes::ELEM_SFNCON4:
                         case ElementTypes::ELEM_SFNCON6:
                         case ElementTypes::ELEM_SFNCON6H:
                         case ElementTypes::ELEM_SFNCON8:
                         case ElementTypes::ELEM_SFNCON8R:
                         case ElementTypes::ELEM_SFNCON10: {
                              const auto iter_constr = m_etype.seek("constraint");

                              if (iter_constr != m_etype.end()) {
                                   ov_constr = m_etype.contents(iter_constr);
                              }
                         } break;
                         default:
                              break;
                         };

                         switch (s) {
                         case STAGE_RESERVE:
                              for (octave_idx_type j = 0; j < ov_fn.numel(); ++j) {
                                   octave_idx_type icurrconstr = 0;

                                   switch (oElemType.type) {
                                   case ElementTypes::ELEM_JOINT:
                                   case ElementTypes::ELEM_THERM_CONSTR:
                                   case ElementTypes::ELEM_ACOUSTIC_CONSTR:
                                        icurrconstr = ov_fn(j).rows();
                                        edof[CS_JOINT].elem_count++;
                                        break;
                                   case ElementTypes::ELEM_SFNCON4:
                                   case ElementTypes::ELEM_SFNCON6:
                                   case ElementTypes::ELEM_SFNCON6H:
                                   case ElementTypes::ELEM_SFNCON8:
                                   case ElementTypes::ELEM_SFNCON8R:
                                   case ElementTypes::ELEM_SFNCON10:
#if HAVE_NLOPT == 1
                                        icurrconstr = SurfToNodeConstrBase::iGetNumDof(ov_constr, j, eDomain);
                                        edof[CS_JOINT].elem_count += ov_fn(j).numel();
#else
                                        error(SurfToNodeConstrBase::szErrCompileWithNlopt);
                                        return retval;
#endif
                                        break;
                                   default:
                                        FEM_ASSERT(false);
                                   }

                                   edof[CS_JOINT].maxdof = std::max(edof[CS_JOINT].maxdof, icurrconstr);
                              }
                              break;
                         case STAGE_FILL:
                              FEM_ASSERT(edof[CS_JOINT].dof.rows() == edof[CS_JOINT].elem_count ||
                                         edof[CS_JOINT].dof.rows() == 0);

                              edof[CS_JOINT].dof.resize(dim_vector(edof[CS_JOINT].elem_count,
                                                                   edof[CS_JOINT].maxdof),
                                                        0);

                              for (octave_idx_type j = 0; j < ov_fn.numel(); ++j) {
                                   octave_idx_type icurrconstr = 0, icurrjoints = 0;

                                   switch (oElemType.type) {
                                   case ElementTypes::ELEM_JOINT:
                                   case ElementTypes::ELEM_THERM_CONSTR:
                                   case ElementTypes::ELEM_ACOUSTIC_CONSTR:
                                        icurrconstr = ov_fn(j).rows();
                                        icurrjoints = 1;
                                        break;
                                   case ElementTypes::ELEM_SFNCON4:
                                   case ElementTypes::ELEM_SFNCON6:
                                   case ElementTypes::ELEM_SFNCON6H:
                                   case ElementTypes::ELEM_SFNCON8:
                                   case ElementTypes::ELEM_SFNCON8R:
                                   case ElementTypes::ELEM_SFNCON10:
#if HAVE_NLOPT == 1
                                        icurrconstr = SurfToNodeConstrBase::iGetNumDof(ov_constr, j, eDomain);
                                        icurrjoints += ov_fn(j).numel();
#else
                                        error(SurfToNodeConstrBase::szErrCompileWithNlopt);
                                        return retval;
#endif
                                        break;
                                   default:
                                        FEM_ASSERT(false);
                                   }

                                   for (octave_idx_type k = 0; k < icurrjoints; ++k) {
                                        for (octave_idx_type l = 0; l < icurrconstr; ++l) {
                                             edof[CS_JOINT].dof.xelem(edof[CS_JOINT].elem_idx, l) = ++icurrdof;
                                             ++inumlambda;
                                        }
                                        edof[CS_JOINT].elem_idx++;
                                   }
                              }

                              m_edof.assign("joints", edof[CS_JOINT].dof);
                         }
                    } break;
                    case ElementTypes::ELEM_RBE3:
                         switch (s) {
                         case STAGE_RESERVE:
                              break;
                         case STAGE_FILL:
                              edof[CS_RBE3].elem_count = m_etype.numel();
                              edof[CS_RBE3].maxdof = 6;

                              edof[CS_RBE3].dof.resize(dim_vector(edof[CS_RBE3].elem_count,
                                                                  edof[CS_RBE3].maxdof),
                                                       0);

                              for (octave_idx_type j = 0; j < edof[CS_RBE3].elem_count; ++j) {
                                   for (octave_idx_type l = 0; l < edof[CS_RBE3].maxdof; ++l) {
                                        edof[CS_RBE3].dof.xelem(j + edof[CS_RBE3].elem_count * l) = ++icurrdof;
                                        ++inumlambda;
                                   }
                              }

                              m_edof.assign("rbe3", edof[CS_RBE3].dof);
                         } break;
                    default:
                         continue;
                    }
               }
          }

          int32NDArray idx_lambda(dim_vector(inumlambda, 1), -1);
          octave_idx_type icurrlambda = 0;

          for (octave_idx_type k = 0; k < CS_COUNT; ++k) {
               const octave_idx_type edofkrows = edof[k].dof.rows();
               const octave_idx_type edofkcols = edof[k].dof.columns();
               for (octave_idx_type i = 0; i < edofkrows; ++i) {
                    for (octave_idx_type j = 0; j < edofkcols; ++j) {
                         const octave_idx_type idxedof = edof[k].dof.xelem(i + edofkrows * j).value();

                         if (idxedof > 0) {
                              idx_lambda.xelem(icurrlambda++) = idxedof;
                         }
                    }
               }
          }

          idx_lambda.sort();

          FEM_ASSERT(icurrlambda == inumlambda);

          if (m_edof.nfields()) {
               dof_map.assign("edof", m_edof);
          }

          if (inumlambda) {
               dof_map.assign("idx_lambda", idx_lambda);
          }

          dof_map.assign("totdof", octave_int32(icurrdof));

          retval.append(dof_map);

     } catch (const std::exception& err) {
          error("%s", err.what());
     }

     return retval;
}

DEFUN_DLD(fem_pre_mesh_constr_surf_to_node, args, nargout,
          "-*- texinfo -*-\n"
          "@deftypefn {} @var{joints} = fem_pre_mesh_constr_surf_to_node(@var{nodes}, @var{elements}, @var{domain})\n"
          "Convert elements of type sfncon4 or sfncon6 to joints.\n"
          "Slave nodes outside @var{elements}.sfncon4.maxdist or @var{elements}.sfncon6.maxdist are ignored.\n\n"
          "@var{nodes} @dots{} The node position matrix of the unconstrained finite element mesh.\n\n"
          "@var{elements}.sfncon4 @dots{} Struct array of elements which are connecting slave nodes to quadrilateral elements.\n\n"
          "@var{elements}.sfncon4.slave @dots{} Array of slave node numbers.\n\n"
          "@var{elements}.sfncon4.master @dots{} Master mesh defined by four node quadrilateral surface elements.\n\n"
          "@var{elements}.sfncon4.maxdist @dots{} Maximum distance between slave nodes and master mesh.\n\n"
          "@var{elements}.sfncon6 @dots{} Struct array of elements which are connecting slave nodes to triangular elements.\n\n"
          "@var{elements}.sfncon6.slave @dots{} Array of slave node numbers.\n\n"
          "@var{elements}.sfncon6.master @dots{} Master mesh defined by six node triangular surface elements.\n\n"
          "@var{elements}.sfncon6.maxdist @dots{} Maximum distance between slave nodes and master mesh.\n\n"
          "@var{joints} @dots{} Struct array of joint elements, one for each slave node within a distance of\n"
          "<@var{maxdist}> from the master mesh.\n\n"
          "@seealso{fem_tests}\n\n"
          "@end deftypefn\n")
{
     octave_value_list retval;
#if HAVE_NLOPT == 1
     const octave_idx_type nargin = args.length();

     if (!(nargin >= 2 && nargin <= 3)) {
          print_usage();
          return retval;
     }

     const Matrix nodes(args(0).matrix_value());

#if OCTAVE_MAJOR_VERSION < 6
     if (error_state) {
          return retval;
     }
#endif

     const octave_scalar_map elements(args(1).scalar_map_value());

#if OCTAVE_MAJOR_VERSION < 6
     if (error_state) {
          return retval;
     }
#endif

     DofMap::DomainType eDomain = DofMap::DO_STRUCTURAL;

     if (nargin > 2) {
          eDomain = static_cast<DofMap::DomainType>(args(2).int_value());
     }

     const octave_idx_type iNodeMaxDofIndex = DofMap::iGetNodeMaxDofIndex(eDomain);

     if (iNodeMaxDofIndex <= 0) {
          error("invalid value for domain=%d", static_cast<int>(eDomain));
          return retval;
     }

     try {
          if (nodes.columns() != 6) {
               throw std::runtime_error("fem_pre_mesh_constr_surf_to_node: invalid number of columns for matrix nodes");
          }

          array<int32NDArray, DofMap::ELEM_TYPE_COUNT> edof;
          array<octave_idx_type, DofMap::ELEM_TYPE_COUNT> dofelemid = {0};

          vector<std::unique_ptr<ElementBlockBase> > rgElemBlocks;

          rgElemBlocks.reserve(2);

          for (octave_idx_type k = ElementTypes::ELEM_SFNCON4; k <= ElementTypes::ELEM_SFNCON10; ++k) {
               constexpr unsigned uFlags = SurfToNodeConstrBase::CF_IGNORE_NODES_OUT_OF_RANGE;
               const ElementTypes::TypeInfo& oElemType = ElementTypes::GetType(k);

               FEM_ASSERT(oElemType.type == k);

               SurfToNodeConstrBase::BuildJoints(nodes, elements, edof, dofelemid, oElemType, uFlags, rgElemBlocks, eDomain);
          }

          octave_idx_type iNumElem = 0;

          for (const auto& oElemBlk: rgElemBlocks) {
               iNumElem += oElemBlk->iGetNumElem();
          }

          constexpr size_t nFields = 2;
          static constexpr const char* const rgFields[nFields] = {"C", "nodes"};
          const string_vector strFields(rgFields, nFields);

          octave_map sElem(dim_vector(1, iNumElem), strFields);
          octave_idx_type idx = 0;

          for (const auto& oElemBlk: rgElemBlocks) {
               oElemBlk->Extract(idx, sElem);
          }

          sElem.resize(dim_vector(1, idx));

          retval.append(sElem);
     } catch (const std::exception& err) {
          error("%s", err.what());
          return retval;
     }
#else
     error(SurfToNodeConstrBase::szErrCompileWithNlopt);
#endif
     return retval;
}

DEFUN_DLD(fem_ass_matrix, args, nargout,
          "-*- texinfo -*-\n"
          "@deftypefn {} [@var{varargout}] = fem_ass_matrix(@var{mesh}, @var{dof_map}, @var{matrix_type})\n"
          "@deftypefnx {} [@dots{}] = fem_ass_matrix(@var{mesh}, @var{dof_map}, @var{matrix_type}, @var{load_case})\n"
          "@deftypefnx {} [@dots{}] = fem_ass_matrix(@var{mesh}, @var{dof_map}, @var{matrix_type}, @var{load_case}, @var{sol})\n"
          "@deftypefnx {} [@dots{}, @var{mat_info}] = fem_ass_matrix(@dots{})\n"
          "@deftypefnx {} [@dots{}, @var{mat_info}, @var{mesh_info}] = fem_ass_matrix(@dots{})\n"
          "This function is the core of the finite element toolkit.\n\n"
          "Assemble all global finite element matrices requested in the array @var{matrix_type} and return it in @var{varargout}.\n\n"
          "@var{mesh} @dots{} Finite element mesh data structure returned from fem_pre_mesh_struct_create or fem_pre_mesh_unstruct_create.\n\n"
          "@var{dof_map} @dots{} Degree of freedom mapping returned from fem_ass_dof_map.\n\n"
          "@var{matrix_type} @dots{} Array of integer numbers identifying the matrix type (e.g. FEM_MAT_*, FEM_VEC_*, FEM_SCA_*).\n\n"
          "@var{load_case} @dots{} Struct Array of loads.\n\n"
          "@var{sol} @dots{} Finite element solution returned from fem_sol_static, fem_sol_modal or fem_post_cms_expand.\n\n"
          "@seealso{fem_tests}\n\n"
          "@end deftypefn\n")
{
     octave_value_list retval;

     const octave_idx_type nargin = args.length();

     if (nargin < 3 || nargin > 5)
     {
          print_usage();
          return retval;
     }

     try {
          const octave_scalar_map mesh(args(0).scalar_map_value());

#if OCTAVE_MAJOR_VERSION < 6
          if (error_state) {
               throw std::runtime_error("fem_ass_matrix: argument mesh must be a scalar struct");
          }
#endif

          const auto it_nodes = mesh.seek("nodes");

          if (it_nodes == mesh.end()) {
               throw std::runtime_error("fem_ass_matrix: missing field mesh.nodes in argument mesh");
          }

          const Matrix nodes(mesh.contents(it_nodes).matrix_value());

#if OCTAVE_MAJOR_VERSION < 6
          if (error_state) {
               throw std::runtime_error("fem_ass_matrix: mesh.nodes must be a real matrix in argument mesh");
          }
#endif
          const octave_idx_type iNumNodes = nodes.rows();

          const auto it_elements = mesh.seek("elements");

          if (it_elements == mesh.end()) {
               throw std::runtime_error("fem_ass_matrix: missing field mesh.elements in argument mesh");
          }

          const octave_scalar_map elements(mesh.contents(it_elements).scalar_map_value());

#if OCTAVE_MAJOR_VERSION < 6
          if (error_state) {
               throw std::runtime_error("fem_ass_matrix: mesh.elements must be a scalar struct in argument mesh");
          }
#endif
          const auto it_materials = mesh.seek("materials");

          if (it_materials == mesh.end()) {
               throw std::runtime_error("fem_ass_matrix: missing field mesh.materials in argument mesh");
          }

          const octave_scalar_map materials(mesh.contents(it_materials).scalar_map_value());

#if OCTAVE_MAJOR_VERSION < 6
          if (error_state) {
               throw std::runtime_error("fem_ass_matrix: mesh.materials must be a scalar struct in argument mesh");
          }
#endif

          const auto it_material_data = mesh.seek("material_data");

          if (it_material_data == mesh.end()) {
               throw std::runtime_error("fem_ass_matrix: missing field mesh.material_data in argument mesh");
          }

          const octave_map material_data(mesh.contents(it_material_data).map_value());

#if OCTAVE_MAJOR_VERSION < 6
          if (error_state) {
               throw std::runtime_error("fem_ass_matrix: mesh.material_data must be a struct array in argument mesh");
          }
#endif

          const octave_scalar_map dof_map(args(1).scalar_map_value());

#if OCTAVE_MAJOR_VERSION < 6
          if (error_state) {
               throw std::runtime_error("fem_ass_matrix: argument dof_map must be a scalar struct");
          }
#endif

          const int32NDArray matrix_type(args(2).int32_array_value());

#if OCTAVE_MAJOR_VERSION < 6
          if (error_state) {
               throw std::runtime_error("fem_ass_matrix: argument matrix_type must be an array of integers");
          }
#endif

          const octave_map load_case(nargin > 3 ? args(3).map_value() : octave_map());

          const octave_idx_type iNumLoads = load_case.numel();

#if OCTAVE_MAJOR_VERSION < 6
          if (error_state) {
               throw std::runtime_error("fem_ass_matrix: argument load case must be a struct array");
          }
#endif

          Element3D::ElementData oElemData(load_case, nodes, elements);

          const octave_scalar_map sol(nargin > 4 ? args(4).scalar_map_value() : octave_scalar_map());

#if OCTAVE_MAJOR_VERSION < 6
          if (error_state) {
               throw std::runtime_error("fem_ass_matrix: argument sol must be scalar struct");
          }
#endif

          if (nodes.columns() != 6) {
               throw std::runtime_error("fem_ass_matrix: invalid number of columns for matrix mesh.nodes in argument mesh");
          }

          const auto iter_ndof = dof_map.seek("ndof");

          if (iter_ndof == dof_map.end()) {
               throw std::runtime_error("fem_ass_matrix: field \"ndof\" not found in argument dof_map");
          }

          const auto iter_domain = dof_map.seek("domain");

          if (iter_domain == dof_map.end()) {
               throw std::runtime_error("fem_ass_matrix: field \"domain\" not found in argument dof_map");
          }

          const auto eDomain = static_cast<DofMap::DomainType>(dof_map.contents(iter_domain).int_value());

          const int32NDArray ndof(dof_map.contents(iter_ndof).int32_array_value());

#if OCTAVE_MAJOR_VERSION < 6
          if (error_state) {
               throw std::runtime_error("fem_ass_matrix: field dof_map.ndof must be an integer array in argument dof_map");
          }
#endif

          const auto iter_totdof = dof_map.seek("totdof");

          if (iter_totdof == dof_map.end()) {
               throw std::runtime_error("fem_ass_matrix: field \"totdof\" not found in argument dof_map");
          }

          const octave_idx_type inumdof = dof_map.contents(iter_totdof).int32_scalar_value();

#if OCTAVE_MAJOR_VERSION < 6
          if (error_state) {
               throw std::runtime_error("fem_ass_matrix: field dof_map.totdof must be an scalar integer in argument dof_map");
          }
#endif

          if (ndof.rows() != iNumNodes || ndof.columns() != DofMap::iGetNodeMaxDofIndex(eDomain)) {
               throw std::runtime_error("fem_ass_matrix: shape of dof_map.ndof is not valid");
          }

          for (octave_idx_type j = 0; j < ndof.columns(); ++j) {
               for (octave_idx_type i = 0; i < iNumNodes; ++i) {
                    octave_idx_type idof = ndof(i, j).value();
                    if (idof > inumdof) {
                         throw std::runtime_error("fem_ass_matrix: invalid index in matrix dof_map.ndof in argument dof_map");
                    }
               }
          }

          const auto iter_edof = dof_map.seek("edof");

          array<int32NDArray, DofMap::ELEM_TYPE_COUNT> edof;
          array<octave_idx_type, DofMap::ELEM_TYPE_COUNT> dofelemid = {0};

          if (iter_edof != dof_map.end()) {
               const octave_scalar_map s_edof(dof_map.contents(iter_edof).scalar_map_value());

#if OCTAVE_MAJOR_VERSION < 6
               if (error_state) {
                    throw std::runtime_error("fem_ass_matrix: dof_map.edof must be a scalar struct in argument dof_map");
               }
#endif

               static constexpr struct DofEntries {
                    DofMap::ElementType type;
                    char name[7];
                    octave_idx_type col_min, col_max;
               } rgDofEntries[] = {
                    {DofMap::ELEM_RBE3, "rbe3", 6, 6},
                    {DofMap::ELEM_JOINT, "joints", 1, -1}
               };

               for (auto k = std::begin(rgDofEntries); k != std::end(rgDofEntries); ++k) {
                    const auto iter_dof = s_edof.seek(k->name);

                    if (iter_dof != s_edof.end()) {
                         const int32NDArray a_edof(s_edof.contents(iter_dof).int32_array_value());

#if OCTAVE_MAJOR_VERSION < 6
                         if (error_state) {
                              std::ostringstream os;
                              os << "fem_ass_matrix: dof_map.edof." << k->name << " must be an integer array in argument dof_map";
                              throw std::runtime_error(os.str());
                         }
#endif

                         if (a_edof.columns() < k->col_min || (k->col_max > 0 && a_edof.columns() > k->col_max)) {
                              std::ostringstream os;
                              os << "fem_ass_matrix: dof_map.edof." << k->name << " numer of columns = " << a_edof.columns() << " not in range [" << k->col_min << ":" << k->col_max << "] in argument dof_map";
                              throw std::runtime_error(os.str());
                         }

                         edof[k->type] = a_edof;
                    }
               }
          }

          for (const auto& edofk: edof) {
               const octave_idx_type edofkrows = edofk.rows();
               const octave_idx_type edofkcols = edofk.columns();

               for (octave_idx_type j = 0; j < edofkcols; ++j) {
                    for (octave_idx_type i = 0; i < edofkrows; ++i) {
                         const octave_idx_type edofkij = edofk.xelem(i + edofkrows * j).value();

                         if (edofkij > inumdof) {
                              throw std::runtime_error("fem_ass_matrix: dof_map.edof dof index out of range in argument dof_map");
                         }
                    }
               }
          }

          const ParallelOptions oParaOpt{dof_map};

          const vector<Material> rgMaterials = Material::ExtractMaterialData(material_data, eDomain);

          DofMap oDof(eDomain, ndof, edof, inumdof);

          PostProcData oSolution{oDof, sol};

          array<bool, ElementTypes::iGetNumTypes()> rgElemUse;

          std::fill(std::begin(rgElemUse), std::end(rgElemUse), false);

          for (octave_idx_type i = 0; i < matrix_type.numel(); ++i) {
               const auto eMatType = static_cast<Element::FemMatrixType>(matrix_type(i).value());

               if (!(eMatType & eDomain)) {
                    throw std::runtime_error("fem_ass_matrix: matrix type is not valid for selected dof_map.domain");
               }

               switch (oDof.GetDomain()) {
               case DofMap::DO_STRUCTURAL:
                    switch (eMatType) {
                    case Element::MAT_STIFFNESS:
                    case Element::MAT_STIFFNESS_SYM:
                    case Element::MAT_STIFFNESS_SYM_L:
                    case Element::MAT_STIFFNESS_TAU0: // FIXME: Needed only for solid elements; should be placed somewhere else
                    case Element::MAT_STIFFNESS_OMEGA:
                    case Element::MAT_STIFFNESS_OMEGA_DOT:
                    case Element::MAT_DAMPING_OMEGA:
                         rgElemUse[ElementTypes::ELEM_RBE3] = true;
                         rgElemUse[ElementTypes::ELEM_JOINT] = true;
                         rgElemUse[ElementTypes::ELEM_SPRING] = true;
                         rgElemUse[ElementTypes::ELEM_SFNCON4] = true;
                         rgElemUse[ElementTypes::ELEM_SFNCON6] = true;
                         rgElemUse[ElementTypes::ELEM_SFNCON6H] = true;
                         rgElemUse[ElementTypes::ELEM_SFNCON8] = true;
                         rgElemUse[ElementTypes::ELEM_SFNCON8R] = true;
                         rgElemUse[ElementTypes::ELEM_SFNCON10] = true;
                         rgElemUse[ElementTypes::ELEM_SFNCON4S] = true;
                         rgElemUse[ElementTypes::ELEM_SFNCON6S] = true;
                         rgElemUse[ElementTypes::ELEM_SFNCON6HS] = true;
                         rgElemUse[ElementTypes::ELEM_SFNCON8S] = true;
                         rgElemUse[ElementTypes::ELEM_SFNCON8RS] = true;
                         rgElemUse[ElementTypes::ELEM_SFNCON10S] = true;
                         [[fallthrough]];

                    case Element::MAT_MASS:
                    case Element::MAT_MASS_SYM:
                    case Element::MAT_MASS_SYM_L:
                    case Element::MAT_MASS_LUMPED:
                    case Element::SCA_TOT_MASS:
                    case Element::VEC_INERTIA_M1:
                    case Element::MAT_INERTIA_J:
                    case Element::MAT_INERTIA_INV3:
                    case Element::MAT_INERTIA_INV4:
                    case Element::MAT_INERTIA_INV5:
                    case Element::MAT_INERTIA_INV8:
                    case Element::MAT_INERTIA_INV9:
                         rgElemUse[ElementTypes::ELEM_RIGID_BODY] = true;
                         [[fallthrough]];

                    case Element::MAT_DAMPING:
                    case Element::MAT_DAMPING_SYM:
                    case Element::MAT_DAMPING_SYM_L:
                         rgElemUse[ElementTypes::ELEM_DASHPOT] = true;
                         [[fallthrough]];

                    case Element::MAT_STIFFNESS_IM:
                         rgElemUse[ElementTypes::ELEM_BEAM2] = true;
                         rgElemUse[ElementTypes::ELEM_SPRING] = true;
                         rgElemUse[ElementTypes::ELEM_SFNCON4S] = true;
                         rgElemUse[ElementTypes::ELEM_SFNCON6S] = true;
                         rgElemUse[ElementTypes::ELEM_SFNCON6HS] = true;
                         rgElemUse[ElementTypes::ELEM_SFNCON8S] = true;
                         rgElemUse[ElementTypes::ELEM_SFNCON8RS] = true;
                         rgElemUse[ElementTypes::ELEM_SFNCON10S] = true;
                         [[fallthrough]];

                    case Element::VEC_STRESS_CAUCH:
                    case Element::VEC_STRAIN_TOTAL:
                    case Element::SCA_STRESS_VMIS:
                    case Element::VEC_COLL_MASS:
                    case Element::VEC_COLL_STIFFNESS:
                         rgElemUse[ElementTypes::ELEM_ISO8] = true;
                         rgElemUse[ElementTypes::ELEM_ISO20] = true;
                         rgElemUse[ElementTypes::ELEM_ISO20R] = true;
                         rgElemUse[ElementTypes::ELEM_ISO27] = true;
                         rgElemUse[ElementTypes::ELEM_PENTA15] = true;
                         rgElemUse[ElementTypes::ELEM_TET10H] = true;
                         rgElemUse[ElementTypes::ELEM_TET10] = true;
                         rgElemUse[ElementTypes::ELEM_TET20] = true;
                         break;

                    case Element::VEC_LOAD_CONSISTENT:
                    case Element::VEC_LOAD_LUMPED:
                         if (iNumLoads == 0) {
                              throw std::runtime_error("fem_ass_matrix: missing argument load_case for matrix_type == FEM_VEC_LOAD_*");
                         }

                         rgElemUse[ElementTypes::ELEM_BEAM2] = true;
                         rgElemUse[ElementTypes::ELEM_RIGID_BODY] = true;
                         rgElemUse[ElementTypes::ELEM_PRESSURE_ISO4] = true;
                         rgElemUse[ElementTypes::ELEM_PRESSURE_TRIA6] = true;
                         rgElemUse[ElementTypes::ELEM_PRESSURE_TRIA6H] = true;
                         rgElemUse[ElementTypes::ELEM_PRESSURE_QUAD8] = true;
                         rgElemUse[ElementTypes::ELEM_PRESSURE_QUAD8R] = true;
                         rgElemUse[ElementTypes::ELEM_PRESSURE_QUAD9] = true;
                         rgElemUse[ElementTypes::ELEM_PRESSURE_TRIA10] = true;
                         rgElemUse[ElementTypes::ELEM_STRUCT_FORCE] = true;

                         rgElemUse[ElementTypes::ELEM_JOINT] = true;
                         // Needed for thermal stress only
                         rgElemUse[ElementTypes::ELEM_ISO8] = true;
                         rgElemUse[ElementTypes::ELEM_ISO20] = true;
                         rgElemUse[ElementTypes::ELEM_ISO20R] = true;
                         rgElemUse[ElementTypes::ELEM_ISO27] = true;
                         rgElemUse[ElementTypes::ELEM_PENTA15] = true;
                         rgElemUse[ElementTypes::ELEM_TET10H] = true;
                         rgElemUse[ElementTypes::ELEM_TET10] = true;
                         rgElemUse[ElementTypes::ELEM_TET20] = true;
                         break;
                    case Element::VEC_SURFACE_AREA:
                         rgElemUse[ElementTypes::ELEM_PRESSURE_ISO4] = true;
                         rgElemUse[ElementTypes::ELEM_PRESSURE_TRIA6] = true;
                         rgElemUse[ElementTypes::ELEM_PRESSURE_TRIA6H] = true;
                         rgElemUse[ElementTypes::ELEM_PRESSURE_QUAD8] = true;
                         rgElemUse[ElementTypes::ELEM_PRESSURE_QUAD8R] = true;
                         rgElemUse[ElementTypes::ELEM_PRESSURE_QUAD9] = true;
                         rgElemUse[ElementTypes::ELEM_PRESSURE_TRIA10] = true;
                         break;
                    default:
                         throw std::runtime_error("fem_ass_matrix: invalid value for argument matrix_type");
                    }
                    break;

               case DofMap::DO_THERMAL:
                    switch (eMatType) {
                    case Element::MAT_THERMAL_COND:
                         rgElemUse[ElementTypes::ELEM_THERM_CONSTR] = true;
                         rgElemUse[ElementTypes::ELEM_SFNCON4] = true;
                         rgElemUse[ElementTypes::ELEM_SFNCON6] = true;
                         rgElemUse[ElementTypes::ELEM_SFNCON6H] = true;
                         rgElemUse[ElementTypes::ELEM_SFNCON8] = true;
                         rgElemUse[ElementTypes::ELEM_SFNCON8R] = true;
                         rgElemUse[ElementTypes::ELEM_SFNCON10] = true;
                         [[fallthrough]];

                    case Element::MAT_HEAT_CAPACITY:
                    case Element::VEC_COLL_HEAT_CAPACITY:
                    case Element::VEC_COLL_THERMAL_COND:
                         rgElemUse[ElementTypes::ELEM_ISO8] = true;
                         rgElemUse[ElementTypes::ELEM_ISO20] = true;
                         rgElemUse[ElementTypes::ELEM_ISO20R] = true;
                         rgElemUse[ElementTypes::ELEM_ISO27] = true;
                         rgElemUse[ElementTypes::ELEM_PENTA15] = true;
                         rgElemUse[ElementTypes::ELEM_TET10H] = true;
                         rgElemUse[ElementTypes::ELEM_TET10] = true;
                         rgElemUse[ElementTypes::ELEM_TET20] = true;
                         break;

                    case Element::VEC_LOAD_THERMAL:
                         if (iNumLoads == 0) {
                              throw std::runtime_error("fem_ass_matrix: missing argument load_case for matrix_type == FEM_VEC_LOAD_THERMAL");
                         }

                         rgElemUse[ElementTypes::ELEM_THERM_CONSTR] = true;
                         rgElemUse[ElementTypes::ELEM_THERM_CONV_ISO4] = true;
                         rgElemUse[ElementTypes::ELEM_THERM_CONV_QUAD8] = true;
                         rgElemUse[ElementTypes::ELEM_THERM_CONV_QUAD8R] = true;
                         rgElemUse[ElementTypes::ELEM_THERM_CONV_QUAD9] = true;
                         rgElemUse[ElementTypes::ELEM_THERM_CONV_TRIA6] = true;
                         rgElemUse[ElementTypes::ELEM_THERM_CONV_TRIA6H] = true;
                         rgElemUse[ElementTypes::ELEM_THERM_CONV_TRIA10] = true;

                         rgElemUse[ElementTypes::ELEM_HEAT_SOURCE_ISO4] = true;
                         rgElemUse[ElementTypes::ELEM_HEAT_SOURCE_QUAD8] = true;
                         rgElemUse[ElementTypes::ELEM_HEAT_SOURCE_QUAD8R] = true;
                         rgElemUse[ElementTypes::ELEM_HEAT_SOURCE_QUAD9] = true;
                         rgElemUse[ElementTypes::ELEM_HEAT_SOURCE_TRIA6] = true;
                         rgElemUse[ElementTypes::ELEM_HEAT_SOURCE_TRIA6H] = true;
                         rgElemUse[ElementTypes::ELEM_HEAT_SOURCE_TRIA10] = true;
                         break;

                    default:
                         throw std::runtime_error("fem_ass_matrix: invalid value for argument matrix_type");
                    }
                    break;

               case DofMap::DO_ACOUSTICS:
                    switch (eMatType) {
                    case Element::VEC_PARTICLE_VELOCITY:
                    case Element::VEC_PARTICLE_VELOCITY_C:
                    case Element::SCA_ACOUSTIC_INTENSITY:
                    case Element::SCA_ACOUSTIC_INTENSITY_C:
                         rgElemUse[ElementTypes::ELEM_ACOUSTIC_BND_ISO4] = true;
                         rgElemUse[ElementTypes::ELEM_ACOUSTIC_BND_QUAD8] = true;
                         rgElemUse[ElementTypes::ELEM_ACOUSTIC_BND_QUAD8R] = true;
                         rgElemUse[ElementTypes::ELEM_ACOUSTIC_BND_QUAD9] = true;
                         rgElemUse[ElementTypes::ELEM_ACOUSTIC_BND_TRIA6] = true;
                         rgElemUse[ElementTypes::ELEM_ACOUSTIC_BND_TRIA6H] = true;
                         rgElemUse[ElementTypes::ELEM_ACOUSTIC_BND_TRIA10] = true;
                         [[fallthrough]];

                    case Element::MAT_MASS_ACOUSTICS_RE:
                    case Element::MAT_MASS_ACOUSTICS_IM:
                    case Element::MAT_STIFFNESS_ACOUSTICS_RE:
                    case Element::MAT_STIFFNESS_ACOUSTICS_IM:
                    case Element::VEC_COLL_MASS_ACOUSTICS:
                    case Element::VEC_COLL_STIFF_ACOUSTICS:
                         rgElemUse[ElementTypes::ELEM_ISO8] = true;
                         rgElemUse[ElementTypes::ELEM_ISO20] = true;
                         rgElemUse[ElementTypes::ELEM_ISO20R] = true;
                         rgElemUse[ElementTypes::ELEM_ISO27] = true;
                         rgElemUse[ElementTypes::ELEM_PENTA15] = true;
                         rgElemUse[ElementTypes::ELEM_TET10H] = true;
                         rgElemUse[ElementTypes::ELEM_TET10] = true;
                         rgElemUse[ElementTypes::ELEM_TET20] = true;
                         break;

                    case Element::VEC_SURFACE_NORMAL_VECTOR:
                         rgElemUse[ElementTypes::ELEM_PARTICLE_VEL_ISO4] = true;
                         rgElemUse[ElementTypes::ELEM_PARTICLE_VEL_QUAD8] = true;
                         rgElemUse[ElementTypes::ELEM_PARTICLE_VEL_QUAD8R] = true;
                         rgElemUse[ElementTypes::ELEM_PARTICLE_VEL_QUAD9] = true;
                         rgElemUse[ElementTypes::ELEM_PARTICLE_VEL_TRIA6] = true;
                         rgElemUse[ElementTypes::ELEM_PARTICLE_VEL_TRIA6H] = true;
                         rgElemUse[ElementTypes::ELEM_PARTICLE_VEL_TRIA10] = true;
                         break;

                    case Element::VEC_LOAD_ACOUSTICS:
                         if (iNumLoads == 0) {
                              throw std::runtime_error("fem_ass_matrix: missing argument load_case for matrix_type == FEM_VEC_LOAD_ACOUSTICS");
                         }

                         rgElemUse[ElementTypes::ELEM_PARTICLE_VEL_ISO4] = true;
                         rgElemUse[ElementTypes::ELEM_PARTICLE_VEL_QUAD8] = true;
                         rgElemUse[ElementTypes::ELEM_PARTICLE_VEL_QUAD8R] = true;
                         rgElemUse[ElementTypes::ELEM_PARTICLE_VEL_QUAD9] = true;
                         rgElemUse[ElementTypes::ELEM_PARTICLE_VEL_TRIA6] = true;
                         rgElemUse[ElementTypes::ELEM_PARTICLE_VEL_TRIA6H] = true;
                         rgElemUse[ElementTypes::ELEM_PARTICLE_VEL_TRIA10] = true;
                         rgElemUse[ElementTypes::ELEM_ACOUSTIC_CONSTR] = true;
                         break;

                    case Element::MAT_DAMPING_ACOUSTICS_RE:
                         rgElemUse[ElementTypes::ELEM_ISO8] = true;
                         rgElemUse[ElementTypes::ELEM_ISO20] = true;
                         rgElemUse[ElementTypes::ELEM_ISO20R] = true;
                         rgElemUse[ElementTypes::ELEM_ISO27] = true;
                         rgElemUse[ElementTypes::ELEM_PENTA15] = true;
                         rgElemUse[ElementTypes::ELEM_TET10H] = true;
                         rgElemUse[ElementTypes::ELEM_TET10] = true;
                         rgElemUse[ElementTypes::ELEM_TET20] = true;
                         rgElemUse[ElementTypes::ELEM_ACOUSTIC_CONSTR] = true;
                         rgElemUse[ElementTypes::ELEM_SFNCON4] = true;
                         rgElemUse[ElementTypes::ELEM_SFNCON6] = true;
                         rgElemUse[ElementTypes::ELEM_SFNCON6H] = true;
                         rgElemUse[ElementTypes::ELEM_SFNCON8] = true;
                         rgElemUse[ElementTypes::ELEM_SFNCON8R] = true;
                         rgElemUse[ElementTypes::ELEM_SFNCON10] = true;
                         [[fallthrough]];

                    case Element::MAT_DAMPING_ACOUSTICS_IM:
                         rgElemUse[ElementTypes::ELEM_ACOUSTIC_IMPE_ISO4] = true;
                         rgElemUse[ElementTypes::ELEM_ACOUSTIC_IMPE_QUAD8] = true;
                         rgElemUse[ElementTypes::ELEM_ACOUSTIC_IMPE_QUAD8R] = true;
                         rgElemUse[ElementTypes::ELEM_ACOUSTIC_IMPE_QUAD9] = true;
                         rgElemUse[ElementTypes::ELEM_ACOUSTIC_IMPE_TRIA6] = true;
                         rgElemUse[ElementTypes::ELEM_ACOUSTIC_IMPE_TRIA6H] = true;
                         rgElemUse[ElementTypes::ELEM_ACOUSTIC_IMPE_TRIA10] = true;
                         break;

                    default:
                         break;
                    }
                    break;

               case DofMap::DO_FLUID_STRUCT:
                    switch (eMatType) {
                    case Element::MAT_STIFFNESS_FLUID_STRUCT_RE:
                         rgElemUse[ElementTypes::ELEM_RBE3] = true;
                         rgElemUse[ElementTypes::ELEM_JOINT] = true;
                         rgElemUse[ElementTypes::ELEM_SPRING] = true;
                         // FIXME: Add support for ELEM_SFNCON*
                         [[fallthrough]];
                    case Element::MAT_MASS_FLUID_STRUCT_RE:
                    case Element::MAT_MASS_FLUID_STRUCT_IM:
                         rgElemUse[ElementTypes::ELEM_BEAM2] = true;
                         rgElemUse[ElementTypes::ELEM_RIGID_BODY] = true;
                         [[fallthrough]];
                    case Element::MAT_STIFFNESS_FLUID_STRUCT_IM:
                    case Element::VEC_COLL_MASS_FLUID_STRUCT:
                    case Element::VEC_COLL_STIFF_FLUID_STRUCT:
                         rgElemUse[ElementTypes::ELEM_ISO8] = true;
                         rgElemUse[ElementTypes::ELEM_ISO20] = true;
                         rgElemUse[ElementTypes::ELEM_ISO20R] = true;
                         rgElemUse[ElementTypes::ELEM_ISO27] = true;
                         rgElemUse[ElementTypes::ELEM_PENTA15] = true;
                         rgElemUse[ElementTypes::ELEM_TET10H] = true;
                         rgElemUse[ElementTypes::ELEM_TET10] = true;
                         rgElemUse[ElementTypes::ELEM_TET20] = true;
                         rgElemUse[ElementTypes::ELEM_SPRING] = true;
                         break;

                    case Element::VEC_LOAD_FLUID_STRUCT:
                         if (iNumLoads == 0) {
                              throw std::runtime_error("fem_ass_matrix: missing argument load_case for matrix_type == FEM_VEC_LOAD_FLUID_STRUCT");
                         }

                         rgElemUse[ElementTypes::ELEM_PRESSURE_ISO4] = true;
                         rgElemUse[ElementTypes::ELEM_PRESSURE_TRIA6] = true;
                         rgElemUse[ElementTypes::ELEM_PRESSURE_TRIA6H] = true;
                         rgElemUse[ElementTypes::ELEM_PRESSURE_TRIA10] = true;
                         rgElemUse[ElementTypes::ELEM_PRESSURE_QUAD8] = true;
                         rgElemUse[ElementTypes::ELEM_PRESSURE_QUAD8R] = true;
                         rgElemUse[ElementTypes::ELEM_PRESSURE_QUAD9] = true;
                         rgElemUse[ElementTypes::ELEM_STRUCT_FORCE] = true;
                         rgElemUse[ElementTypes::ELEM_JOINT] = true;
                         // Needed for thermal stress only
                         rgElemUse[ElementTypes::ELEM_ISO8] = true;
                         rgElemUse[ElementTypes::ELEM_ISO20] = true;
                         rgElemUse[ElementTypes::ELEM_ISO20R] = true;
                         rgElemUse[ElementTypes::ELEM_ISO27] = true;
                         rgElemUse[ElementTypes::ELEM_PENTA15] = true;
                         rgElemUse[ElementTypes::ELEM_TET10H] = true;
                         rgElemUse[ElementTypes::ELEM_TET10] = true;
                         rgElemUse[ElementTypes::ELEM_TET20] = true;

                         rgElemUse[ElementTypes::ELEM_PARTICLE_VEL_ISO4] = true;
                         rgElemUse[ElementTypes::ELEM_PARTICLE_VEL_QUAD8] = true;
                         rgElemUse[ElementTypes::ELEM_PARTICLE_VEL_QUAD8R] = true;
                         rgElemUse[ElementTypes::ELEM_PARTICLE_VEL_QUAD9] = true;
                         rgElemUse[ElementTypes::ELEM_PARTICLE_VEL_TRIA6] = true;
                         rgElemUse[ElementTypes::ELEM_PARTICLE_VEL_TRIA6H] = true;
                         rgElemUse[ElementTypes::ELEM_PARTICLE_VEL_TRIA10] = true;
                         rgElemUse[ElementTypes::ELEM_ACOUSTIC_CONSTR] = true;
                         break;

                    case Element::VEC_SURFACE_AREA:
                         rgElemUse[ElementTypes::ELEM_PRESSURE_ISO4] = true;
                         rgElemUse[ElementTypes::ELEM_PRESSURE_TRIA6] = true;
                         rgElemUse[ElementTypes::ELEM_PRESSURE_TRIA6H] = true;
                         rgElemUse[ElementTypes::ELEM_PRESSURE_QUAD8] = true;
                         rgElemUse[ElementTypes::ELEM_PRESSURE_QUAD8R] = true;
                         rgElemUse[ElementTypes::ELEM_PRESSURE_QUAD9] = true;
                         rgElemUse[ElementTypes::ELEM_PRESSURE_TRIA10] = true;
                         break;

                    case Element::MAT_DAMPING_FLUID_STRUCT_RE:
                         rgElemUse[ElementTypes::ELEM_ISO8] = true;
                         rgElemUse[ElementTypes::ELEM_ISO20] = true;
                         rgElemUse[ElementTypes::ELEM_ISO20R] = true;
                         rgElemUse[ElementTypes::ELEM_ISO27] = true;
                         rgElemUse[ElementTypes::ELEM_PENTA15] = true;
                         rgElemUse[ElementTypes::ELEM_TET10H] = true;
                         rgElemUse[ElementTypes::ELEM_TET10] = true;
                         rgElemUse[ElementTypes::ELEM_TET20] = true;
                         rgElemUse[ElementTypes::ELEM_FLUID_STRUCT_ISO4] = true;
                         rgElemUse[ElementTypes::ELEM_FLUID_STRUCT_QUAD8] = true;
                         rgElemUse[ElementTypes::ELEM_FLUID_STRUCT_QUAD8R] = true;
                         rgElemUse[ElementTypes::ELEM_FLUID_STRUCT_QUAD9] = true;
                         rgElemUse[ElementTypes::ELEM_FLUID_STRUCT_TRIA6] = true;
                         rgElemUse[ElementTypes::ELEM_FLUID_STRUCT_TRIA6H] = true;
                         rgElemUse[ElementTypes::ELEM_FLUID_STRUCT_TRIA10] = true;
                         rgElemUse[ElementTypes::ELEM_ACOUSTIC_CONSTR] = true;
                         rgElemUse[ElementTypes::ELEM_DASHPOT] = true;
                         // FIXME: Add support for ELEM_SFNCON*
                         [[fallthrough]];

                    case Element::MAT_DAMPING_FLUID_STRUCT_IM:
                         rgElemUse[ElementTypes::ELEM_ACOUSTIC_IMPE_ISO4] = true;
                         rgElemUse[ElementTypes::ELEM_ACOUSTIC_IMPE_QUAD8] = true;
                         rgElemUse[ElementTypes::ELEM_ACOUSTIC_IMPE_QUAD8R] = true;
                         rgElemUse[ElementTypes::ELEM_ACOUSTIC_IMPE_QUAD9] = true;
                         rgElemUse[ElementTypes::ELEM_ACOUSTIC_IMPE_TRIA6] = true;
                         rgElemUse[ElementTypes::ELEM_ACOUSTIC_IMPE_TRIA6H] = true;
                         rgElemUse[ElementTypes::ELEM_ACOUSTIC_IMPE_TRIA10] = true;
                         break;

                    case Element::VEC_SURFACE_NORMAL_VECTOR:
                         rgElemUse[ElementTypes::ELEM_PARTICLE_VEL_ISO4] = true;
                         rgElemUse[ElementTypes::ELEM_PARTICLE_VEL_QUAD8] = true;
                         rgElemUse[ElementTypes::ELEM_PARTICLE_VEL_QUAD8R] = true;
                         rgElemUse[ElementTypes::ELEM_PARTICLE_VEL_QUAD9] = true;
                         rgElemUse[ElementTypes::ELEM_PARTICLE_VEL_TRIA6] = true;
                         rgElemUse[ElementTypes::ELEM_PARTICLE_VEL_TRIA6H] = true;
                         rgElemUse[ElementTypes::ELEM_PARTICLE_VEL_TRIA10] = true;
                         break;

                    case Element::VEC_PARTICLE_VELOCITY:
                    case Element::VEC_PARTICLE_VELOCITY_C:
                    case Element::SCA_ACOUSTIC_INTENSITY:
                    case Element::SCA_ACOUSTIC_INTENSITY_C:
                         rgElemUse[ElementTypes::ELEM_ISO8] = true;
                         rgElemUse[ElementTypes::ELEM_ISO20] = true;
                         rgElemUse[ElementTypes::ELEM_ISO20R] = true;
                         rgElemUse[ElementTypes::ELEM_ISO27] = true;
                         rgElemUse[ElementTypes::ELEM_PENTA15] = true;
                         rgElemUse[ElementTypes::ELEM_TET10H] = true;
                         rgElemUse[ElementTypes::ELEM_TET10] = true;
                         rgElemUse[ElementTypes::ELEM_TET20] = true;
                         rgElemUse[ElementTypes::ELEM_ACOUSTIC_BND_ISO4] = true;
                         rgElemUse[ElementTypes::ELEM_ACOUSTIC_BND_QUAD8] = true;
                         rgElemUse[ElementTypes::ELEM_ACOUSTIC_BND_QUAD8R] = true;
                         rgElemUse[ElementTypes::ELEM_ACOUSTIC_BND_QUAD9] = true;
                         rgElemUse[ElementTypes::ELEM_ACOUSTIC_BND_TRIA6] = true;
                         rgElemUse[ElementTypes::ELEM_ACOUSTIC_BND_TRIA6H] = true;
                         rgElemUse[ElementTypes::ELEM_ACOUSTIC_BND_TRIA10] = true;
                         break;

                    default:
                         break;
                    }
                    break;

               default:
                    throw std::runtime_error("fem_ass_matrix: invalid value for dof_map.domain");
               }
          }

          vector<std::unique_ptr<ElementBlockBase> > rgElemBlocks;

          rgElemBlocks.reserve(ElementTypes::iGetNumTypes());

          for (octave_idx_type k = 0; k < ElementTypes::iGetNumTypes(); ++k) {
               const ElementTypes::TypeInfo& oElemType = ElementTypes::GetType(k);

               if (!rgElemUse[oElemType.type]) {
                    continue;
               }

               switch (oElemType.type) {
               case ElementTypes::ELEM_ISO8:
               case ElementTypes::ELEM_ISO20:
               case ElementTypes::ELEM_ISO20R:
               case ElementTypes::ELEM_ISO27:
               case ElementTypes::ELEM_PENTA15:
               case ElementTypes::ELEM_TET10H:
               case ElementTypes::ELEM_TET10:
               case ElementTypes::ELEM_TET20: {
                    const auto iter_elem = elements.seek(oElemType.name);

                    if (iter_elem == elements.end()) {
                         continue;
                    }

                    const int32NDArray elem_nodes = elements.contents(iter_elem).int32_array_value();

#if OCTAVE_MAJOR_VERSION < 6
                    if (error_state) {
                         throw std::runtime_error("fem_ass_matrix: mesh.elements."s + oElemType.name
                                                  + " must be an array of integers in argument mesh");
                    }
#endif
                    const octave_idx_type iNumElem = elem_nodes.rows();
                    const octave_idx_type iNumNodesElem = elem_nodes.columns();

                    if (iNumNodesElem < oElemType.min_nodes) {
                         throw std::runtime_error("invalid number of nodes for element type "s + oElemType.name
                                                  +  " in argument mesh");
                    }

                    if (oElemType.max_nodes > 0 && iNumNodesElem > oElemType.max_nodes) {
                         throw std::runtime_error("invalid number of nodes for element type "s + oElemType.name
                                                  + " in argument mesh");
                    }

                    for (octave_idx_type j = 0; j < iNumNodesElem; ++j) {
                         for (octave_idx_type i = 0; i < iNumElem; ++i) {
                              const octave_idx_type inode = elem_nodes.xelem(i + iNumElem * j);
                              if (inode < 1 || inode > iNumNodes) {
                                   throw std::runtime_error("invalid node index for element type "s
                                                            + oElemType.name + " in argument mesh");
                              }
                         }
                    }

                    const auto iter_mat = materials.seek(oElemType.name);

                    if (iter_mat == materials.end()) {
                         throw std::runtime_error("fem_ass_matrix: material not defined for all element types in argument mesh");
                    }

                    const int32NDArray elem_mat = materials.contents(iter_mat).int32_array_value();

                    if (elem_mat.numel() != iNumElem) {
                         throw std::runtime_error("fem_ass_matrix: invalid size for matrix mesh.materials."s + oElemType.name + " in argument mesh");
                    }

                    for (octave_idx_type i = 0; i < iNumElem; ++i) {
                         octave_idx_type imaterial = elem_mat.xelem(i);

                         if (imaterial <= 0 || imaterial > material_data.numel()) {
                              throw std::runtime_error("fem_ass_matrix: invalid index in matrix mesh.materials in argument mesh");
                         }
                    }

                    const auto iter_PMLelem = oElemData.oPML.oData.seek(oElemType.name);

                    if (iter_PMLelem != oElemData.oPML.oData.end()) {
                         const octave_scalar_map ov_PMLelem = oElemData.oPML.oData.contents(iter_PMLelem).scalar_map_value();

                         const auto iter_f = ov_PMLelem.seek("f");

                         if (iter_f == ov_PMLelem.end()) {
                              throw std::runtime_error("fem_ass_matrix: missing field mesh.perfectly_matched_layers."s + oElemType.name + ".f");
                         }

                         const octave_value ov_fel = ov_PMLelem.contents(iter_f);

                         if (!((ov_fel.isreal() || ov_fel.iscomplex()) && ov_fel.is_matrix_type())) {
                              throw std::runtime_error("fem_ass_matrix: mesh.perfectly_matched_layers."s + oElemType.name + " must be a real matrix");
                         }

                         if (ov_fel.rows() != 3 || ov_fel.columns() != iNumNodesElem || ov_fel.numel() != 3 * elem_nodes.numel()) {
                              throw std::runtime_error("fem_ass_matrix: invalid size for matrix mesh.perfectly_matched_layers."s + oElemType.name + " in argument mesh");
                         }

                         oElemData.oPML.rgElem[oElemType.type].f = ov_fel.complex_array_value();

                         const auto iter_e1 = ov_PMLelem.seek("e1");

                         if (iter_e1 == ov_PMLelem.end()) {
                              throw std::runtime_error("fem_ass_matrix: missing field mesh.perfectly_matched_layers."s + oElemType.name + ".e1");
                         }

                         const octave_value ov_e1 = ov_PMLelem.contents(iter_e1);

                         if (!(ov_e1.isreal() && ov_e1.is_matrix_type())) {
                              throw std::runtime_error("fem_ass_matrix: mesh.perfectly_matched_layers."s + oElemType.name + ".e1 must be a real matrix");
                         }

                         if (ov_e1.rows() != 3 || ov_e1.columns() != iNumNodesElem || ov_e1.numel() != 3 * elem_nodes.numel()) {
                              throw std::runtime_error("fem_ass_matrix: invalid size for matrix mesh.perfectly_matched_layers."s + oElemType.name + ".e1 in argument mesh");
                         }

                         oElemData.oPML.rgElem[oElemType.type].e1 = ov_e1.array_value();

                         const auto iter_e2 = ov_PMLelem.seek("e2");

                         if (iter_e2 == ov_PMLelem.end()) {
                              throw std::runtime_error("fem_ass_matrix: missing field mesh.perfectly_matched_layers."s + oElemType.name + ".e2");
                         }

                         const octave_value ov_e2 = ov_PMLelem.contents(iter_e2);

                         if (!(ov_e2.isreal() && ov_e2.is_matrix_type())) {
                              throw std::runtime_error("fem_ass_matrix: mesh.perfectly_matched_layers."s + oElemType.name + ".e2 must be a real matrix");
                         }

                         if (ov_e2.rows() != 3 || ov_e2.columns() != iNumNodesElem || ov_e2.numel() != 3 * elem_nodes.numel()) {
                              throw std::runtime_error("fem_ass_matrix: invalid size for matrix mesh.perfectly_matched_layers."s + oElemType.name + ".e2 in argument mesh");
                         }

                         oElemData.oPML.rgElem[oElemType.type].e2 = ov_e2.array_value();
                    }

                    switch (oElemType.type) {
                    case ElementTypes::ELEM_ISO8:
                         rgElemBlocks.emplace_back(new ElementBlock<Iso8>(oElemType.type, elem_nodes, nodes, 3, elem_mat, rgMaterials, oElemData));
                         break;

                    case ElementTypes::ELEM_ISO20:
                         rgElemBlocks.emplace_back(new ElementBlock<Iso20>(oElemType.type, elem_nodes, nodes, 3, elem_mat, rgMaterials, oElemData));
                         break;

                    case ElementTypes::ELEM_ISO20R:
                         rgElemBlocks.emplace_back(new ElementBlock<Iso20r>(oElemType.type, elem_nodes, nodes, 3, elem_mat, rgMaterials, oElemData));
                         break;

                    case ElementTypes::ELEM_ISO27:
                         rgElemBlocks.emplace_back(new ElementBlock<Iso27>(oElemType.type, elem_nodes, nodes, 3, elem_mat, rgMaterials, oElemData));
                         break;

                    case ElementTypes::ELEM_PENTA15:
                         rgElemBlocks.emplace_back(new ElementBlock<Penta15>(oElemType.type, elem_nodes, nodes, 3, elem_mat, rgMaterials, oElemData));
                         break;

                    case ElementTypes::ELEM_TET10H:
                         rgElemBlocks.emplace_back(new ElementBlock<Tet10h>(oElemType.type, elem_nodes, nodes, 3, elem_mat, rgMaterials, oElemData));
                         break;

                    case ElementTypes::ELEM_TET10:
                         rgElemBlocks.emplace_back(new ElementBlock<Tet10>(oElemType.type, elem_nodes, nodes, 3, elem_mat, rgMaterials, oElemData));
                         break;

                    case ElementTypes::ELEM_TET20:
                         rgElemBlocks.emplace_back(new ElementBlock<Tet20>(oElemType.type, elem_nodes, nodes, 3, elem_mat, rgMaterials, oElemData));
                         break;

                    default:
                         throw std::logic_error("fem_ass_matrix: invalid element type");
                    }
               } break;
               case ElementTypes::ELEM_BEAM2:
               case ElementTypes::ELEM_RIGID_BODY:
               case ElementTypes::ELEM_RBE3:
               case ElementTypes::ELEM_JOINT:
               case ElementTypes::ELEM_SPRING:
               case ElementTypes::ELEM_DASHPOT:
               case ElementTypes::ELEM_THERM_CONSTR:
               case ElementTypes::ELEM_ACOUSTIC_CONSTR: {
                    const auto iter_elem = elements.seek(oElemType.name);

                    if (iter_elem == elements.end()) {
                         continue;
                    }

                    const octave_map s_elem(elements.contents(iter_elem).map_value());

#if OCTAVE_MAJOR_VERSION < 6
                    if (error_state) {
                         throw std::runtime_error("fem_ass_matrix: mesh.elements."s + oElemType.name + " must be an struct array in argument mesh");
                    }
#endif
                    const auto iter_nodes = s_elem.seek("nodes");

                    if (iter_nodes == s_elem.end()) {
                         throw std::runtime_error("fem_ass_matrix: missing field mesh.elements."s + oElemType.name + ".nodes in argument mesh");
                    }

                    const Cell ov_nodes(s_elem.contents(iter_nodes));

                    FEM_ASSERT(ov_nodes.numel() == s_elem.numel());

                    int32NDArray elem_mat;
                    Cell ov_section, ov_e2;

                    if (oElemType.type == ElementTypes::ELEM_BEAM2) {
                         const auto iter_mat = materials.seek("beam2");

                         if (iter_mat == materials.end()) {
                              throw std::runtime_error("fem_ass_matrix: missing field mesh.materials.beam2 in argument mesh");
                         }

                         const octave_value ov_elem_mat = materials.contents(iter_mat);

                         if (!(ov_elem_mat.columns() == 1 && ov_elem_mat.isinteger())) {
                              throw std::runtime_error("fem_ass_matrix: mesh.materials.beam2 must be an integer column vector");
                         }

                         elem_mat = ov_elem_mat.int32_array_value();

                         if (elem_mat.numel() != s_elem.numel()) {
                              throw std::runtime_error("fem_ass_matrix: numel(mesh.materials.beam2) does not match numel(mesh.elements.beam2)");
                         }

                         const auto iter_sect = s_elem.seek("section");

                         if (iter_sect == s_elem.end()) {
                              throw std::runtime_error("fem_ass_matrix: missing field mesh.elements.beam2.section in argument mesh");
                         }

                         ov_section = s_elem.contents(iter_sect);

                         FEM_ASSERT(ov_section.numel() == s_elem.numel());

                         const auto iter_e2 = s_elem.seek("e2");

                         if (iter_e2 == s_elem.end()) {
                              throw std::runtime_error("fem_ass_matrix: missing field mesh.elements.beam2.e2 in argument mesh");
                         }

                         ov_e2 = s_elem.contents(iter_e2);

                         FEM_ASSERT(ov_e2.numel() == s_elem.numel());
                    }

                    Cell cell_m, cell_J, cell_lcg;

                    if (oElemType.type == ElementTypes::ELEM_RIGID_BODY) {
                         const auto iter_m = s_elem.seek("m");

                         if (iter_m == s_elem.end()) {
                              throw std::runtime_error("fem_ass_matrix: missing field mesh.elements."s + oElemType.name + ".m in argument mesh");
                         }

                         cell_m = s_elem.contents(iter_m);

                         const auto iter_J = s_elem.seek("J");

                         if (iter_J == s_elem.end()) {
                              throw std::runtime_error("fem_ass_matrix: missing field mesh.elements."s + oElemType.name + ".J in argument mesh");
                         }

                         cell_J = s_elem.contents(iter_J);

                         const auto iter_lcg = s_elem.seek("lcg");

                         if (iter_lcg == s_elem.end()) {
                              throw std::runtime_error("fem_ass_matrix: missing field mesh.elements."s + oElemType.name + ".lcg in argument mesh");
                         }

                         cell_lcg = s_elem.contents(iter_lcg);
                    }

                    Cell ov_A, ov_Scale;

                    if (oElemType.type == ElementTypes::ELEM_JOINT ||
                        oElemType.type == ElementTypes::ELEM_THERM_CONSTR ||
                        oElemType.type == ElementTypes::ELEM_ACOUSTIC_CONSTR) {
                         const auto iter_C = s_elem.seek("C");

                         if (iter_C == s_elem.end()) {
                              throw std::runtime_error("fem_ass_matrix: missing field mesh.elements."s + oElemType.name + ".C in argument mesh");
                         }

                         ov_A = s_elem.contents(iter_C);

                         FEM_ASSERT(ov_A.numel() == s_elem.numel());

                         const auto iter_Scale = s_elem.seek("scale");

                         if (iter_Scale != s_elem.end()) {
                              ov_Scale = s_elem.contents(iter_Scale);
                         }
                    }

                    bool bRealA = true;

                    if (oElemType.type == ElementTypes::ELEM_SPRING) {
                         const auto iter_K = s_elem.seek("K");

                         if (iter_K == s_elem.end()) {
                              throw std::runtime_error("fem_ass_matrix: missing field mesh.elements."s + oElemType.name + ".K in argument mesh");
                         }

                         ov_A = s_elem.contents(iter_K);

                         FEM_ASSERT(ov_A.numel() == s_elem.numel());
                    }

                    if (oElemType.type == ElementTypes::ELEM_DASHPOT) {
                         const auto iter_D = s_elem.seek("D");

                         if (iter_D == s_elem.end()) {
                              throw std::runtime_error("fem_ass_matrix: missing field mesh.elements."s + oElemType.name + ".D in argument mesh");
                         }

                         ov_A = s_elem.contents(iter_D);

                         FEM_ASSERT(ov_A.numel() == s_elem.numel());
                    }

                    Cell ov_weight;

                    if (oElemType.type == ElementTypes::ELEM_RBE3) {
                         const auto iter_weight = s_elem.seek("weight");

                         if (iter_weight != s_elem.end()) {
                              ov_weight = s_elem.contents(iter_weight);

                              FEM_ASSERT(ov_weight.numel() == s_elem.numel());
                         }
                    }

                    std::unique_ptr<ElementBlockBase> pElem;

                    switch (oElemType.type) {
                    case ElementTypes::ELEM_BEAM2:
                         pElem.reset(new ElementBlock<ElemBeam2>(oElemType.type, s_elem.numel()));
                         break;
                    case ElementTypes::ELEM_RIGID_BODY:
                         pElem.reset(new ElementBlock<ElemRigidBody>(oElemType.type, s_elem.numel()));
                         break;
                    case ElementTypes::ELEM_RBE3:
                         pElem.reset(new ElementBlock<ElemRBE3>(oElemType.type, s_elem.numel()));
                         break;
                    case ElementTypes::ELEM_JOINT:
                    case ElementTypes::ELEM_THERM_CONSTR:
                    case ElementTypes::ELEM_ACOUSTIC_CONSTR:
                         pElem.reset(new ElementBlock<ElemJoint>(oElemType.type, s_elem.numel()));
                         break;
                    case ElementTypes::ELEM_SPRING:
                    case ElementTypes::ELEM_DASHPOT: {
                         for (octave_idx_type i = 0; i < ov_A.numel(); ++i) {
                              if (!(ov_A.xelem(i).isreal() && ov_A.xelem(i).isreal())) {
                                   bRealA = false;
                                   break;
                              }
                         }

                         switch (oElemType.type) {
                         case ElementTypes::ELEM_SPRING:
                              if (bRealA) {
                                   pElem.reset(new ElementBlock<ElemSpring<double, Matrix>>(oElemType.type, s_elem.numel()));
                              } else {
                                   pElem.reset(new ElementBlock<ElemSpring<std::complex<double>, ComplexMatrix>>(oElemType.type, s_elem.numel()));
                              }
                              break;
                         case ElementTypes::ELEM_DASHPOT:
                              if (bRealA) {
                                   pElem.reset(new ElementBlock<ElemDashpot<double, Matrix>>(oElemType.type, s_elem.numel()));
                              } else {
                                   pElem.reset(new ElementBlock<ElemDashpot<std::complex<double>, ComplexMatrix>>(oElemType.type, s_elem.numel()));
                              }
                              break;
                         default:
                              FEM_ASSERT(false);
                         }
                    } break;
                    default:
                         FEM_ASSERT(false);
                    }

                    for (octave_idx_type i = 0; i < s_elem.numel(); ++i) {
                         const int32NDArray elem_nodes(ov_nodes(i).int32_array_value());

#if OCTAVE_MAJOR_VERSION < 6
                         if (error_state) {
                              throw std::runtime_error("fem_ass_matrix: mesh.elements"s + oElemType.name + ".nodes must be an array of integers in argument mesh");
                         }
#endif

                         if (elem_nodes.columns() < oElemType.min_nodes) {
                              throw std::runtime_error("fem_ass_matrix: invalid number of nodes for element type "s + oElemType.name + " in argument mesh");
                         }

                         if (oElemType.max_nodes > 0 && elem_nodes.columns() > oElemType.max_nodes) {
                              throw std::runtime_error("fem_ass_matrix: invalid number of nodes for element type "s + oElemType.name + " in argument mesh");
                         }

                         if (elem_nodes.rows() != 1) {
                              throw std::runtime_error("fem_ass_matrix: invalid number of rows in node matrix for element type "s + oElemType.name + " in argument mesh");
                         }

                         for (octave_idx_type j = 0; j < elem_nodes.columns(); ++j) {
                              octave_idx_type inode = elem_nodes.xelem(j);

                              if (inode < 1 || inode > iNumNodes) {
                                   throw std::runtime_error("fem_ass_matrix: invalid node index for element type "s + oElemType.name + " in argument mesh");
                              }
                         }

                         Matrix X(6, elem_nodes.columns());

                         for (octave_idx_type j = 0; j < elem_nodes.columns(); ++j) {
                              for (octave_idx_type l = 0; l < 6; ++l) {
                                   X.xelem(l + 6 * j) = nodes.xelem(elem_nodes(j).value() - 1 + iNumNodes * l);
                              }
                         }

                         switch (oElemType.type) {
                         case ElementTypes::ELEM_BEAM2: {
                              const octave_idx_type imaterial(elem_mat.xelem(i).value() - 1);

#if OCTAVE_MAJOR_VERSION < 6
                              if (error_state) {
                                   throw std::runtime_error("fem_ass_matrix: mesh.materials.beam2 must be an integer array");
                              }
#endif

                              if (imaterial < 0 || static_cast<size_t>(imaterial) >= rgMaterials.size()) {
                                   throw std::runtime_error("fem_ass_matrix: mesh.materials.beam2 out of range");
                              }

                              const Material* const material = &rgMaterials[imaterial];

                              const octave_scalar_map m_section(ov_section(i).scalar_map_value());

#if OCTAVE_MAJOR_VERSION < 6
                              if (error_state) {
                                   throw std::runtime_error("fem_ass_matrix: mesh.elements.beam2.section must be a scalar struct");
                              }
#endif

                              auto iterA = m_section.seek("A");

                              if (iterA == m_section.end()) {
                                   throw std::runtime_error("fem_ass_matrix: missing field mesh.elements.beam2.section.A");
                              }

                              auto iterAy = m_section.seek("Ay");

                              if (iterAy == m_section.end()) {
                                   throw std::runtime_error("fem_ass_matrix: missing field mesh.elements.beam2.section.Ay");
                              }

                              auto iterAz = m_section.seek("Az");

                              if (iterAz == m_section.end()) {
                                   throw std::runtime_error("fem_ass_matrix: missing field mesh.elements.beam2.section.Az");
                              }

                              auto iterIt = m_section.seek("It");

                              if (iterIt == m_section.end()) {
                                   throw std::runtime_error("fem_ass_matrix: missing field mesh.elements.beam2.section.It");
                              }

                              auto iterIy = m_section.seek("Iy");

                              if (iterIy == m_section.end()) {
                                   throw std::runtime_error("fem_ass_matrix: missing field mesh.elements.beam2.section.Iy");
                              }

                              auto iterIz = m_section.seek("Iz");

                              if (iterIz == m_section.end()) {
                                   throw std::runtime_error("fem_ass_matrix: missing field mesh.elements.beam2.section.Iz");
                              }

                              BeamCrossSection section;

                              section.A = m_section.contents(iterA).scalar_value();

#if OCTAVE_MAJOR_VERSION < 6
                              if (error_state) {
                                   throw std::runtime_error("fem_ass_matrix: mesh.elements.beam2.section.A must be a real scalar");
                              }
#endif

                              section.Ay = m_section.contents(iterAy).scalar_value();

#if OCTAVE_MAJOR_VERSION < 6
                              if (error_state) {
                                   throw std::runtime_error("fem_ass_matrix: mesh.elements.beam2.section.Ay must be a real scalar");
                              }
#endif

                              section.Az = m_section.contents(iterAz).scalar_value();

#if OCTAVE_MAJOR_VERSION < 6
                              if (error_state) {
                                   throw std::runtime_error("fem_ass_matrix: mesh.elements.beam2.section.Az must be a real scalar");
                              }
#endif

                              section.It = m_section.contents(iterIt).scalar_value();

#if OCTAVE_MAJOR_VERSION < 6
                              if (error_state) {
                                   throw std::runtime_error("fem_ass_matrix: mesh.elements.beam2.section.It must be a real scalar");
                              }
#endif

                              section.Iy = m_section.contents(iterIy).scalar_value();

#if OCTAVE_MAJOR_VERSION < 6
                              if (error_state) {
                                   throw std::runtime_error("fem_ass_matrix: mesh.elements.beam2.section.Iy must be a real scalar");
                              }
#endif

                              section.Iz = m_section.contents(iterIz).scalar_value();

#if OCTAVE_MAJOR_VERSION < 6
                              if (error_state) {
                                   throw std::runtime_error("fem_ass_matrix: mesh.elements.beam2.section.Iz must be a real scalar");
                              }
#endif

                              const ColumnVector e2(ov_e2(i).column_vector_value());

#if OCTAVE_MAJOR_VERSION < 6
                              if (error_state) {
                                   throw std::runtime_error("fem_ass_matrix: mesh.elements.beam2.e2 must be a column vector");
                              }
#endif

                              if (e2.rows() != 3) {
                                   throw std::runtime_error("fem_ass_matrix: mesh.elements.beam2.e2 must be 3x1 matrix");
                              }

                              pElem->Insert<ElemBeam2>(i + 1, X, material, elem_nodes, section, e2, oElemData.oGravity.g);
                         } break;
                         case ElementTypes::ELEM_RIGID_BODY: {
                              const octave_value ov_m = cell_m.xelem(i);

                              if (!ov_m.is_real_scalar()) {
                                   throw std::runtime_error("fem_ass_matrix: mesh.elements.body.m must be a real scalar");
                              }

                              const double m = ov_m.scalar_value();

                              const octave_value ov_J = cell_J.xelem(i);

                              if (!(ov_J.is_matrix_type() && ov_J.isreal() && ov_J.rows() == 3 && ov_J.columns() == 3)) {
                                   throw std::runtime_error("fem_ass_matrix: mesh.elements.body.J must be a symmetric 3x3 matrix");
                              }

                              const Matrix J = ov_J.matrix_value();

                              if (!J.issymmetric()) {
                                   throw std::runtime_error("fem_ass_matrix: mesh.elements.body.J must be a symmetric 3x3 matrix");
                              }

                              const octave_value ov_lcg = cell_lcg.xelem(i);

                              if (!(ov_lcg.is_matrix_type() && ov_lcg.isreal() && ov_lcg.rows() == 3 && ov_lcg.columns() == 1)) {
                                   throw std::runtime_error("fem_ass_matrix: mesh.elements.body.lcg must be a real 3x1 vector");
                              }

                              const ColumnVector lcg = ov_lcg.column_vector_value();

                              pElem->Insert<ElemRigidBody>(i + 1, X, nullptr, elem_nodes, m, J, lcg, oElemData.oGravity.g);
                         } break;
                         case ElementTypes::ELEM_RBE3: {
                              RowVector weight;

                              if (ov_weight.numel() > i) {
                                   weight = ov_weight(i).row_vector_value();

#if OCTAVE_MAJOR_VERSION < 6
                                   if (error_state) {
                                        throw std::runtime_error("fem_ass_matrix: mesh.elements.rbe3.weight must be a row vector in argument mesh");
                                   }
#endif
                              } else {
                                   weight.resize(elem_nodes.numel(), 1.);
                              }

                              if (weight.numel() != elem_nodes.numel() - 1) {
                                   throw std::runtime_error("fem_ass_matrix: numel(mesh.elements.rbe3.weight) does not match numel(mesh.elements.rbe3.nodes) - 1 in argument mesh");
                              }

                              pElem->Insert<ElemRBE3>(++dofelemid[oElemType.dof_type], X, nullptr, elem_nodes, weight);
                         } break;
                         case ElementTypes::ELEM_JOINT:
                         case ElementTypes::ELEM_THERM_CONSTR:
                         case ElementTypes::ELEM_ACOUSTIC_CONSTR: {
                              const Matrix C(ov_A.xelem(i).matrix_value());

#if OCTAVE_MAJOR_VERSION < 6
                              if (error_state) {
                                   throw std::runtime_error("fem_ass_matrix: mesh.elements."s + oElemType.name + ".C must be a real matrix in argument mesh");
                              }
#endif

                              const octave_idx_type iNumConstrEq = C.rows();
                              const octave_idx_type iNumDofConstr = C.columns();

                              const double dScale = ov_Scale.numel() ? ov_Scale.xelem(i).scalar_value() : 1.;

#if OCTAVE_MAJOR_VERSION < 6
                              if (error_state) {
                                   throw std::runtime_error("fem_ass_matrix: mesh.elements."s + oElemType.name + ".scale must be a real scalar in argument mesh");
                              }
#endif

                              const octave_idx_type iNumDofNodeMax = ElemJoint::iGetNumDofNodeMax(oElemType.type);

                              if (iNumConstrEq < 1 || iNumConstrEq > edof[oElemType.dof_type].columns() || iNumDofConstr != iNumDofNodeMax * elem_nodes.columns() || iNumConstrEq > iNumDofConstr) {
                                   throw std::runtime_error("fem_ass_matrix: invalid size for field elements."s + oElemType.name + ".C");
                              }

                              Matrix U(iNumConstrEq, iNumLoads, 0.); // By default displacement is set to zero

                              const auto iter_joints = load_case.seek(oElemType.name);

                              if (iter_joints != load_case.end()) {
                                   const Cell ov_joints = load_case.contents(iter_joints);

                                   for (octave_idx_type m = 0; m < iNumLoads; ++m) {
                                        if (ov_joints(m).isempty()) {
                                             continue;
                                        }

                                        const octave_map s_joints(ov_joints(m).map_value());

#if OCTAVE_MAJOR_VERSION < 6
                                        if (error_state) {
                                             throw std::runtime_error("fem_ass_matrix: load_case."s + oElemType.name + " must be an struct array in argument load_case");
                                        }
#endif

                                        enum ConstrVarName {
                                             CS_VAR_NAME_STRUCT,
                                             CS_VAR_NAME_THERMAL,
                                             CS_VAR_NAME_ACOUSTIC,
                                             CS_VAR_NAME_EMPTY
                                        } idx = CS_VAR_NAME_EMPTY;

                                        static constexpr char var_name[4][6] = {"U", "theta", "p", ""};

                                        switch (oElemType.type) {
                                        case ElementTypes::ELEM_JOINT:
                                             idx = CS_VAR_NAME_STRUCT;
                                             break;
                                        case ElementTypes::ELEM_THERM_CONSTR:
                                             idx = CS_VAR_NAME_THERMAL;
                                             break;
                                        case ElementTypes::ELEM_ACOUSTIC_CONSTR:
                                             idx = CS_VAR_NAME_ACOUSTIC;
                                             break;
                                        default:
                                             FEM_ASSERT(0);
                                        }

                                        const auto iter_U = s_joints.seek(var_name[idx]);

                                        if (iter_U == s_elem.end()) {
                                             throw std::runtime_error("fem_ass_matrix: missing field load_case."s + oElemType.name + ".U in argument load case");
                                        }

                                        const Cell ov_U(s_joints.contents(iter_U));

                                        if (ov_U.numel() != s_elem.numel()) {
                                             throw std::runtime_error("fem_ass_matrix: load_case."s + oElemType.name + " must have the same size like mesh.elements." + oElemType.name + " in argument load case");
                                        }

                                        const ColumnVector Uk(ov_U(i).column_vector_value());

#if OCTAVE_MAJOR_VERSION < 6
                                        if (error_state) {
                                             throw std::runtime_error("fem_ass_matrix: load_case."s + oElemType.name + ".U must be a real column vector");
                                        }
#endif

                                        if (Uk.rows() != iNumConstrEq) {
                                             throw std::runtime_error("fem_ass_matrix: load_case."s + oElemType.name + ".U must be a real column vector of the same number of rows like mesh.elements." + oElemType.name + ".C in argument load_case");
                                        }

                                        for (octave_idx_type l = 0; l < iNumConstrEq; ++l) {
                                             U.xelem(l + iNumConstrEq * m) = Uk.xelem(l);
                                        }
                                   }
                              }

                              pElem->Insert<ElemJoint>(++dofelemid[oElemType.dof_type], X, nullptr, elem_nodes, C, U, oDof.GetDomain(), dScale);
                         } break;
                         case ElementTypes::ELEM_SPRING:
                         case ElementTypes::ELEM_DASHPOT: {
                              if (bRealA) {
                                   const Matrix A(ov_A.xelem(i).matrix_value());

                                   switch (oElemType.type) {
                                   case ElementTypes::ELEM_SPRING:
                                        pElem->Insert<ElemSpring<double, Matrix>>(i + 1, X, nullptr, elem_nodes, A);
                                        break;
                                   case ElementTypes::ELEM_DASHPOT:
                                        pElem->Insert<ElemDashpot<double, Matrix>>(i + 1, X, nullptr, elem_nodes, A);
                                        break;
                                   default:
                                        FEM_ASSERT(false);
                                   }
                              } else {
                                   const ComplexMatrix A(ov_A.xelem(i).complex_matrix_value());

                                   switch (oElemType.type) {
                                   case ElementTypes::ELEM_SPRING:
                                        pElem->Insert<ElemSpring<std::complex<double>, ComplexMatrix>>(i + 1, X, nullptr, elem_nodes, A);
                                        break;
                                   case ElementTypes::ELEM_DASHPOT:
                                        pElem->Insert<ElemDashpot<std::complex<double>, ComplexMatrix>>(i + 1, X, nullptr, elem_nodes, A);
                                        break;
                                   default:
                                        FEM_ASSERT(false);
                                   }
                              }
                         } break;
                         default:
                              FEM_ASSERT(false);
                         }
                    }

                    if (oElemType.dof_type != DofMap::ELEM_NODOF && dofelemid[oElemType.dof_type] > edof[oElemType.dof_type].rows()) {
                         throw std::runtime_error("fem_ass_matrix: dof_map.edof is not consistent with elements in argument dof_map");
                    }

                    rgElemBlocks.emplace_back(std::move(pElem));
               } break;
               case ElementTypes::ELEM_SFNCON4:
               case ElementTypes::ELEM_SFNCON6:
               case ElementTypes::ELEM_SFNCON6H:
               case ElementTypes::ELEM_SFNCON8:
               case ElementTypes::ELEM_SFNCON8R:
               case ElementTypes::ELEM_SFNCON10: {
#if HAVE_NLOPT == 1
                    constexpr unsigned uFlags = SurfToNodeConstrBase::CF_ELEM_DOF_PRE_ALLOCATED;
                    SurfToNodeConstrBase::BuildJoints(nodes, elements, edof, dofelemid, oElemType, uFlags, rgElemBlocks, oDof.GetDomain());
#else
                    throw std::runtime_error(SurfToNodeConstrBase::szErrCompileWithNlopt);
#endif
               } break;
               case ElementTypes::ELEM_SFNCON4S:
               case ElementTypes::ELEM_SFNCON6S:
               case ElementTypes::ELEM_SFNCON6HS:
               case ElementTypes::ELEM_SFNCON8S:
               case ElementTypes::ELEM_SFNCON8RS:
               case ElementTypes::ELEM_SFNCON10S: {
#if HAVE_NLOPT == 1
                    constexpr unsigned uFlags = SurfToNodeConstrBase::CF_ELEM_DOF_PRE_ALLOCATED;
                    SurfToNodeConstrBase::BuildContacts(nodes, elements, oElemType, uFlags, rgElemBlocks, oDof.GetDomain());
#else
                    throw std::runtime_error(SurfToNodeConstrBase::szErrCompileWithNlopt);
#endif
               } break;
               case ElementTypes::ELEM_PRESSURE_ISO4:
                    InsertPressureElem<PressureLoadIso4>(oElemType.type, nodes, load_case, oElemType.name, oElemType.max_nodes, rgElemBlocks);
                    break;
               case ElementTypes::ELEM_PRESSURE_QUAD8:
                    InsertPressureElem<PressureLoadQuad8>(oElemType.type, nodes, load_case, oElemType.name, oElemType.max_nodes, rgElemBlocks);
                    break;
               case ElementTypes::ELEM_PRESSURE_QUAD8R:
                    InsertPressureElem<PressureLoadQuad8r>(oElemType.type, nodes, load_case, oElemType.name, oElemType.max_nodes, rgElemBlocks);
                    break;
               case ElementTypes::ELEM_PRESSURE_QUAD9:
                    InsertPressureElem<PressureLoadQuad9>(oElemType.type, nodes, load_case, oElemType.name, oElemType.max_nodes, rgElemBlocks);
                    break;
               case ElementTypes::ELEM_PRESSURE_TRIA6:
                    InsertPressureElem<PressureLoadTria6>(oElemType.type, nodes, load_case, oElemType.name, oElemType.max_nodes, rgElemBlocks);
                    break;
               case ElementTypes::ELEM_PRESSURE_TRIA6H:
                    InsertPressureElem<PressureLoadTria6H>(oElemType.type, nodes, load_case, oElemType.name, oElemType.max_nodes, rgElemBlocks);
                    break;
               case ElementTypes::ELEM_PRESSURE_TRIA10:
                    InsertPressureElem<PressureLoadTria10>(oElemType.type, nodes, load_case, oElemType.name, oElemType.max_nodes, rgElemBlocks);
                    break;
               case ElementTypes::ELEM_FLUID_STRUCT_ISO4:
                    InsertFluidStructElem<FluidStructInteractIso4>(oElemType.type, nodes, elements, oElemType.name, oElemType.max_nodes, rgElemBlocks);
                    break;
               case ElementTypes::ELEM_FLUID_STRUCT_QUAD8:
                    InsertFluidStructElem<FluidStructInteractQuad8>(oElemType.type, nodes, elements, oElemType.name, oElemType.max_nodes, rgElemBlocks);
                    break;
               case ElementTypes::ELEM_FLUID_STRUCT_QUAD8R:
                    InsertFluidStructElem<FluidStructInteractQuad8r>(oElemType.type, nodes, elements, oElemType.name, oElemType.max_nodes, rgElemBlocks);
                    break;
               case ElementTypes::ELEM_FLUID_STRUCT_QUAD9:
                    InsertFluidStructElem<FluidStructInteractQuad9>(oElemType.type, nodes, elements, oElemType.name, oElemType.max_nodes, rgElemBlocks);
                    break;
               case ElementTypes::ELEM_FLUID_STRUCT_TRIA6:
                    InsertFluidStructElem<FluidStructInteractTria6>(oElemType.type, nodes, elements, oElemType.name, oElemType.max_nodes, rgElemBlocks);
                    break;
               case ElementTypes::ELEM_FLUID_STRUCT_TRIA6H:
                    InsertFluidStructElem<FluidStructInteractTria6H>(oElemType.type, nodes, elements, oElemType.name, oElemType.max_nodes, rgElemBlocks);
                    break;
               case ElementTypes::ELEM_FLUID_STRUCT_TRIA10:
                    InsertFluidStructElem<FluidStructInteractTria10>(oElemType.type, nodes, elements, oElemType.name, oElemType.max_nodes, rgElemBlocks);
                    break;
               case ElementTypes::ELEM_STRUCT_FORCE: {
                    const auto iter_loads = load_case.seek("loads");
                    const auto iter_loaded_nodes = load_case.seek("loaded_nodes");

                    if (iter_loads != load_case.end() && iter_loaded_nodes != load_case.end()) {
                         const Cell cell_loads = load_case.contents(iter_loads);
                         const Cell cell_loaded_nodes = load_case.contents(iter_loaded_nodes);

                         FEM_ASSERT(cell_loads.numel() == cell_loaded_nodes.numel());
                         FEM_ASSERT(cell_loads.numel() == iNumLoads);

                         octave_idx_type iNumForces = 0;

                         for (octave_idx_type i = 0; i < cell_loads.numel(); ++i) {
                              if (cell_loads(i).numel()) {
                                   ++iNumForces;
                              }
                         }

                         if (iNumForces) {
                              std::unique_ptr<ElementBlock<StructForce> > pElem{new ElementBlock<StructForce>(oElemType.type, iNumForces)};

                              for (octave_idx_type i = 0; i < cell_loads.numel(); ++i) {
                                   if (cell_loads(i).numel()) {
                                        const Matrix loads = cell_loads(i).matrix_value();
#if OCTAVE_MAJOR_VERSION < 6
                                        if (error_state) {
                                             throw std::runtime_error("fem_ass_matrix: field load_case.loads must be a real matrix in argument load_case");
                                        }
#endif

                                        if (loads.columns() != 3 && loads.columns() != 6) {
                                             throw std::runtime_error("fem_ass_matrix: load_case.loads must be a n x 3 or n x 6 matrix in argument load_case");
                                        }

                                        const int32NDArray loaded_nodes = cell_loaded_nodes(i).int32_array_value();

#if OCTAVE_MAJOR_VERSION < 6
                                        if (error_state) {
                                             throw std::runtime_error("fem_ass_matrix: field load_case.loaded_nodes must be an integer matrix in argument load_case");
                                        }
#endif

                                        if (loaded_nodes.columns() != 1 || loaded_nodes.rows() != loads.rows()) {
                                             throw std::runtime_error("fem_ass_matrix: load_case.loaded_nodes must be a column vector with the same number of rows like load_case.loads");
                                        }

                                        Matrix X(3, loaded_nodes.rows());

                                        for (octave_idx_type l = 0; l < X.columns(); ++l) {
                                             octave_idx_type inode = loaded_nodes(l).value() - 1;

                                             if (inode < 0 || inode >= nodes.rows()) {
                                                  throw std::runtime_error("fem_ass_matrix: node index out of range in load_case.pressure.elements in argument load_case");
                                             }

                                             for (octave_idx_type m = 0; m < 3; ++m) {
                                                  X.xelem(m + 3 * l) = nodes.xelem(inode, m);
                                             }
                                        }

                                        pElem->Insert(i + 1, X, nullptr, loaded_nodes, i + 1, loads);
                                   }
                              }

                              rgElemBlocks.emplace_back(std::move(pElem));
                         }
                    }
               } break;
               case ElementTypes::ELEM_THERM_CONV_ISO4:
                    InsertThermalConvElem<ThermalConvectionBCIso4>(oElemType.type, nodes, elements, load_case, oElemType.name, oElemType.max_nodes, rgElemBlocks);
                    break;
               case ElementTypes::ELEM_THERM_CONV_QUAD8:
                    InsertThermalConvElem<ThermalConvectionBCQuad8>(oElemType.type, nodes, elements, load_case, oElemType.name, oElemType.max_nodes, rgElemBlocks);
                    break;
               case ElementTypes::ELEM_THERM_CONV_QUAD8R:
                    InsertThermalConvElem<ThermalConvectionBCQuad8r>(oElemType.type, nodes, elements, load_case, oElemType.name, oElemType.max_nodes, rgElemBlocks);
                    break;
               case ElementTypes::ELEM_THERM_CONV_QUAD9:
                    InsertThermalConvElem<ThermalConvectionBCQuad9>(oElemType.type, nodes, elements, load_case, oElemType.name, oElemType.max_nodes, rgElemBlocks);
                    break;
               case ElementTypes::ELEM_THERM_CONV_TRIA6:
                    InsertThermalConvElem<ThermalConvectionBCTria6>(oElemType.type, nodes, elements, load_case, oElemType.name, oElemType.max_nodes, rgElemBlocks);
                    break;
               case ElementTypes::ELEM_THERM_CONV_TRIA6H:
                    InsertThermalConvElem<ThermalConvectionBCTria6H>(oElemType.type, nodes, elements, load_case, oElemType.name, oElemType.max_nodes, rgElemBlocks);
                    break;
               case ElementTypes::ELEM_THERM_CONV_TRIA10:
                    InsertThermalConvElem<ThermalConvectionBCTria10>(oElemType.type, nodes, elements, load_case, oElemType.name, oElemType.max_nodes, rgElemBlocks);
                    break;
               case ElementTypes::ELEM_HEAT_SOURCE_ISO4:
                    InsertHeatSourceElem<HeatSourceIso4>(oElemType.type, nodes, load_case, oElemType.name, oElemType.max_nodes, rgElemBlocks);
                    break;
               case ElementTypes::ELEM_HEAT_SOURCE_QUAD8:
                    InsertHeatSourceElem<HeatSourceQuad8>(oElemType.type, nodes, load_case, oElemType.name, oElemType.max_nodes, rgElemBlocks);
                    break;
               case ElementTypes::ELEM_HEAT_SOURCE_QUAD8R:
                    InsertHeatSourceElem<HeatSourceQuad8r>(oElemType.type, nodes, load_case, oElemType.name, oElemType.max_nodes, rgElemBlocks);
                    break;
               case ElementTypes::ELEM_HEAT_SOURCE_QUAD9:
                    InsertHeatSourceElem<HeatSourceQuad9>(oElemType.type, nodes, load_case, oElemType.name, oElemType.max_nodes, rgElemBlocks);
                    break;
               case ElementTypes::ELEM_HEAT_SOURCE_TRIA6:
                    InsertHeatSourceElem<HeatSourceTria6>(oElemType.type, nodes, load_case, oElemType.name, oElemType.max_nodes, rgElemBlocks);
                    break;
               case ElementTypes::ELEM_HEAT_SOURCE_TRIA6H:
                    InsertHeatSourceElem<HeatSourceTria6H>(oElemType.type, nodes, load_case, oElemType.name, oElemType.max_nodes, rgElemBlocks);
                    break;
               case ElementTypes::ELEM_HEAT_SOURCE_TRIA10:
                    InsertHeatSourceElem<HeatSourceTria10>(oElemType.type, nodes, load_case, oElemType.name, oElemType.max_nodes, rgElemBlocks);
                    break;
               case ElementTypes::ELEM_PARTICLE_VEL_ISO4:
                    InsertParticleVelocityBC<ParticleVelocityBCIso4>(oElemType.type, nodes, elements, rgMaterials, materials, load_case, oElemType.name, oElemType.max_nodes, rgElemBlocks, oDof.GetDomain());
                    break;
               case ElementTypes::ELEM_PARTICLE_VEL_QUAD8:
                    InsertParticleVelocityBC<ParticleVelocityBCQuad8>(oElemType.type, nodes, elements, rgMaterials, materials, load_case, oElemType.name, oElemType.max_nodes, rgElemBlocks, oDof.GetDomain());
                    break;
               case ElementTypes::ELEM_PARTICLE_VEL_QUAD8R:
                    InsertParticleVelocityBC<ParticleVelocityBCQuad8r>(oElemType.type, nodes, elements, rgMaterials, materials, load_case, oElemType.name, oElemType.max_nodes, rgElemBlocks, oDof.GetDomain());
                    break;
               case ElementTypes::ELEM_PARTICLE_VEL_QUAD9:
                    InsertParticleVelocityBC<ParticleVelocityBCQuad9>(oElemType.type, nodes, elements, rgMaterials, materials, load_case, oElemType.name, oElemType.max_nodes, rgElemBlocks, oDof.GetDomain());
                    break;
               case ElementTypes::ELEM_PARTICLE_VEL_TRIA6:
                    InsertParticleVelocityBC<ParticleVelocityBCTria6>(oElemType.type, nodes, elements, rgMaterials, materials, load_case, oElemType.name, oElemType.max_nodes, rgElemBlocks, oDof.GetDomain());
                    break;
               case ElementTypes::ELEM_PARTICLE_VEL_TRIA6H:
                    InsertParticleVelocityBC<ParticleVelocityBCTria6H>(oElemType.type, nodes, elements, rgMaterials, materials, load_case, oElemType.name, oElemType.max_nodes, rgElemBlocks, oDof.GetDomain());
                    break;
               case ElementTypes::ELEM_PARTICLE_VEL_TRIA10:
                    InsertParticleVelocityBC<ParticleVelocityBCTria10>(oElemType.type, nodes, elements, rgMaterials, materials, load_case, oElemType.name, oElemType.max_nodes, rgElemBlocks, oDof.GetDomain());
                    break;
               case ElementTypes::ELEM_ACOUSTIC_IMPE_ISO4:
                    InsertAcousticImpedanceBC<AcousticImpedanceBCIso4>(oElemType.type, nodes, elements, rgMaterials, materials, oElemType.name, oElemType.max_nodes, rgElemBlocks, oDof.GetDomain());
                    break;
               case ElementTypes::ELEM_ACOUSTIC_IMPE_QUAD8:
                    InsertAcousticImpedanceBC<AcousticImpedanceBCQuad8>(oElemType.type, nodes, elements, rgMaterials, materials, oElemType.name, oElemType.max_nodes, rgElemBlocks, oDof.GetDomain());
                    break;
               case ElementTypes::ELEM_ACOUSTIC_IMPE_QUAD8R:
                    InsertAcousticImpedanceBC<AcousticImpedanceBCQuad8r>(oElemType.type, nodes, elements, rgMaterials, materials, oElemType.name, oElemType.max_nodes, rgElemBlocks, oDof.GetDomain());
                    break;
               case ElementTypes::ELEM_ACOUSTIC_IMPE_QUAD9:
                    InsertAcousticImpedanceBC<AcousticImpedanceBCQuad9>(oElemType.type, nodes, elements, rgMaterials, materials, oElemType.name, oElemType.max_nodes, rgElemBlocks, oDof.GetDomain());
                    break;
               case ElementTypes::ELEM_ACOUSTIC_IMPE_TRIA6:
                    InsertAcousticImpedanceBC<AcousticImpedanceBCTria6>(oElemType.type, nodes, elements, rgMaterials, materials, oElemType.name, oElemType.max_nodes, rgElemBlocks, oDof.GetDomain());
                    break;
               case ElementTypes::ELEM_ACOUSTIC_IMPE_TRIA6H:
                    InsertAcousticImpedanceBC<AcousticImpedanceBCTria6H>(oElemType.type, nodes, elements, rgMaterials, materials, oElemType.name, oElemType.max_nodes, rgElemBlocks, oDof.GetDomain());
                    break;
               case ElementTypes::ELEM_ACOUSTIC_IMPE_TRIA10:
                    InsertAcousticImpedanceBC<AcousticImpedanceBCTria10>(oElemType.type, nodes, elements, rgMaterials, materials, oElemType.name, oElemType.max_nodes, rgElemBlocks, oDof.GetDomain());
                    break;
               case ElementTypes::ELEM_ACOUSTIC_BND_ISO4:
                    InsertAcousticBoundary<AcousticBoundaryIso4>(oElemType.type, nodes, elements, rgMaterials, materials, oElemType.name, oElemType.max_nodes, rgElemBlocks);
                    break;
               case ElementTypes::ELEM_ACOUSTIC_BND_QUAD8:
                    InsertAcousticBoundary<AcousticBoundaryQuad8>(oElemType.type, nodes, elements, rgMaterials, materials, oElemType.name, oElemType.max_nodes, rgElemBlocks);
                    break;
               case ElementTypes::ELEM_ACOUSTIC_BND_QUAD8R:
                    InsertAcousticBoundary<AcousticBoundaryQuad8r>(oElemType.type, nodes, elements, rgMaterials, materials, oElemType.name, oElemType.max_nodes, rgElemBlocks);
                    break;
               case ElementTypes::ELEM_ACOUSTIC_BND_QUAD9:
                    InsertAcousticBoundary<AcousticBoundaryQuad9>(oElemType.type, nodes, elements, rgMaterials, materials, oElemType.name, oElemType.max_nodes, rgElemBlocks);
                    break;
               case ElementTypes::ELEM_ACOUSTIC_BND_TRIA6:
                    InsertAcousticBoundary<AcousticBoundaryTria6>(oElemType.type, nodes, elements, rgMaterials, materials, oElemType.name, oElemType.max_nodes, rgElemBlocks);
                    break;
               case ElementTypes::ELEM_ACOUSTIC_BND_TRIA6H:
                    InsertAcousticBoundary<AcousticBoundaryTria6H>(oElemType.type, nodes, elements, rgMaterials, materials, oElemType.name, oElemType.max_nodes, rgElemBlocks);
                    break;
               case ElementTypes::ELEM_ACOUSTIC_BND_TRIA10:
                    InsertAcousticBoundary<AcousticBoundaryTria10>(oElemType.type, nodes, elements, rgMaterials, materials, oElemType.name, oElemType.max_nodes, rgElemBlocks);
                    break;
               default:
                    FEM_ASSERT(false);
               }
          }

          for (octave_idx_type j = 0; j < ElementTypes::iGetNumTypes(); ++j) {
               const auto eDofType = ElementTypes::GetType(j).dof_type;
               if (rgElemUse[j] &&
                   eDofType != DofMap::ELEM_NODOF &&
                   dofelemid[eDofType] != edof[eDofType].rows()) {
                    throw std::runtime_error("fem_ass_matrix: dof_map.edof is not consistent with mesh.elements in argument dof_map");
               }
          }

          octave_idx_type iMaxWorkSpaceSize = 0;

          for (octave_idx_type j = 0; j < matrix_type.numel(); ++j) {
               const Element::FemMatrixType eMatType = static_cast<Element::FemMatrixType>(matrix_type(j).value());
               octave_idx_type iWorkSpaceSize = 0;

               for (const auto& pElemBlock: rgElemBlocks) {
                    iWorkSpaceSize += pElemBlock->iGetWorkSpaceSize(eMatType);
               }

               iMaxWorkSpaceSize = std::max(iMaxWorkSpaceSize, iWorkSpaceSize);
          }

          MatrixAss oMatAss(iMaxWorkSpaceSize);
          MeshInfo oMeshInfo;

          for (octave_idx_type i = 0; i < matrix_type.numel(); ++i) {
               const Element::FemMatrixType eMatType = static_cast<Element::FemMatrixType>(matrix_type.xelem(i).value());

               switch (eMatType) {
               case Element::MAT_STIFFNESS:
               case Element::MAT_STIFFNESS_SYM:
               case Element::MAT_STIFFNESS_SYM_L:
               case Element::MAT_STIFFNESS_IM:
               case Element::MAT_STIFFNESS_TAU0:
               case Element::MAT_STIFFNESS_OMEGA:
               case Element::MAT_STIFFNESS_OMEGA_DOT:
               case Element::MAT_DAMPING_OMEGA:
               case Element::MAT_STIFFNESS_FLUID_STRUCT_RE:
               case Element::MAT_STIFFNESS_FLUID_STRUCT_IM:
               case Element::MAT_MASS:
               case Element::MAT_MASS_SYM:
               case Element::MAT_MASS_SYM_L:
               case Element::MAT_MASS_LUMPED:
               case Element::MAT_MASS_FLUID_STRUCT_RE:
               case Element::MAT_MASS_FLUID_STRUCT_IM:
               case Element::MAT_DAMPING:
               case Element::MAT_DAMPING_SYM:
               case Element::MAT_DAMPING_SYM_L:
               case Element::MAT_DAMPING_FLUID_STRUCT_RE:
               case Element::MAT_DAMPING_FLUID_STRUCT_IM:
               case Element::VEC_LOAD_CONSISTENT:
               case Element::VEC_LOAD_LUMPED:
               case Element::VEC_LOAD_FLUID_STRUCT:
               case Element::MAT_THERMAL_COND:
               case Element::MAT_HEAT_CAPACITY:
               case Element::VEC_LOAD_THERMAL:
               case Element::MAT_STIFFNESS_ACOUSTICS_RE:
               case Element::MAT_STIFFNESS_ACOUSTICS_IM:
               case Element::MAT_MASS_ACOUSTICS_RE:
               case Element::MAT_MASS_ACOUSTICS_IM:
               case Element::MAT_DAMPING_ACOUSTICS_RE:
               case Element::MAT_DAMPING_ACOUSTICS_IM:
               case Element::VEC_LOAD_ACOUSTICS: {
                    oMatAss.Reset(eMatType);

                    bool bMatInfo = false;

                    for (const auto& pElemBlock: rgElemBlocks) {
                         const bool bNeedMatInfo = pElemBlock->bNeedMatrixInfo(eMatType);

                         if (!bMatInfo && bNeedMatInfo) {
                              oMatAss.UpdateMatrixInfo(oDof);
                              bMatInfo = true;
                         }

                         pElemBlock->Assemble(oMatAss, oMeshInfo, oDof, eMatType, oParaOpt);
                    }

                    oMatAss.Finish();

                    retval.append(oMatAss.Assemble(oDof, iNumLoads));

                    if (eMatType & Element::MAT_UPDATE_INFO_ALWAYS) {
                         oMatAss.UpdateMatrixInfo(oDof);
                    }
               } break;
               case Element::SCA_TOT_MASS: {
                    double dMass = 0.;

                    for (const auto& pElemBlock: rgElemBlocks) {
                         dMass += pElemBlock->dGetMass();
                    }

                    retval.append(dMass);
               } break;
               case Element::VEC_INERTIA_M1:
               case Element::MAT_INERTIA_J:
               case Element::MAT_INERTIA_INV3:
               case Element::MAT_INERTIA_INV4:
               case Element::MAT_INERTIA_INV5:
               case Element::MAT_INERTIA_INV8:
               case Element::MAT_INERTIA_INV9: {
                    const octave_idx_type iNumModes = oSolution.GetNumSteps();
                    bool bNeedSolution = false;
                    dim_vector mat_dim;

                    switch (eMatType) {
                    case Element::VEC_INERTIA_M1:
                         mat_dim = dim_vector(3, 1);
                         break;

                    case Element::MAT_INERTIA_J:
                         mat_dim = dim_vector(3, 3);
                         break;

                    case Element::MAT_INERTIA_INV3:
                    case Element::MAT_INERTIA_INV4:
                         mat_dim = dim_vector(3, iNumModes);
                         bNeedSolution = true;
                         break;

                    case Element::MAT_INERTIA_INV5:
                         mat_dim = dim_vector(3, iNumModes, iNumModes);
                         bNeedSolution = true;
                         break;

                    case Element::MAT_INERTIA_INV8:
                         mat_dim = dim_vector(3, 3, iNumModes);
                         bNeedSolution = true;
                         break;

                    case Element::MAT_INERTIA_INV9:
                         mat_dim = dim_vector(3, 3, iNumModes, iNumModes);
                         bNeedSolution = true;
                         break;
                    default:
                         FEM_ASSERT(false);
                    }

                    if (bNeedSolution) {
                         if (nargin <= 4) {
                              throw std::runtime_error("fem_ass_matrix: argument sol is not optional for selected matrix type in argument matrix_type");
                         }
                    }

                    PostProcData::FieldTypeReal eFieldType;

                    switch (eMatType) {
                    case Element::VEC_INERTIA_M1:
                         eFieldType = PostProcData::VEC_G_STRUCT_INERTIA_M1_RE;
                         break;
                    case Element::MAT_INERTIA_J:
                         eFieldType = PostProcData::MAT_G_STRUCT_INERTIA_J_RE;
                         break;
                    case Element::MAT_INERTIA_INV3:
                         eFieldType = PostProcData::MAT_G_STRUCT_INERTIA_INV3_RE;
                         break;
                    case Element::MAT_INERTIA_INV4:
                         eFieldType = PostProcData::MAT_G_STRUCT_INERTIA_INV4_RE;
                         break;
                    case Element::MAT_INERTIA_INV5:
                         eFieldType = PostProcData::MAT_G_STRUCT_INERTIA_INV5_RE;
                         break;
                    case Element::MAT_INERTIA_INV8:
                         eFieldType = PostProcData::MAT_G_STRUCT_INERTIA_INV8_RE;
                         break;
                    case Element::MAT_INERTIA_INV9:
                         eFieldType = PostProcData::MAT_G_STRUCT_INERTIA_INV9_RE;
                         break;
                    default:
                         FEM_ASSERT(0);
                         throw std::logic_error("fem_ass_matrix: unexpected matrix type");
                    }

                    oSolution.SetField(eFieldType, ElementTypes::ELEM_TYPE_UNKNOWN, NDArray(mat_dim, 0.));

                    for (const auto& pElemBlock: rgElemBlocks) {
                         pElemBlock->PostProcElem(eMatType, oSolution, oParaOpt);
                    }

                    NDArray mat = oSolution.GetField(eFieldType, ElementTypes::ELEM_TYPE_UNKNOWN);

                    switch (eMatType) {
                    case Element::MAT_INERTIA_J:
                         FEM_ASSERT(mat.ndims() == 2);
                         FEM_ASSERT(mat.rows() == 3);
                         FEM_ASSERT(mat.columns() == 3);

                         for (octave_idx_type k = 1; k < mat.rows(); ++k) {
                              for (octave_idx_type j = 0; j < k; ++j) {
                                   mat.xelem(k + 3 * j) = mat.xelem(j + 3 * k);
                              }
                         }
                         break;

                    default:
                         break;
                    }

                    retval.append(mat);
               } break;
               case Element::VEC_COLL_MASS:
               case Element::VEC_COLL_STIFFNESS:
               case Element::VEC_COLL_HEAT_CAPACITY:
               case Element::VEC_COLL_THERMAL_COND:
               case Element::VEC_COLL_MASS_ACOUSTICS:
               case Element::VEC_COLL_STIFF_ACOUSTICS:
               case Element::VEC_COLL_MASS_FLUID_STRUCT:
               case Element::VEC_COLL_STIFF_FLUID_STRUCT: {
                    octave_scalar_map mapCollocPoints;

                    for (const auto& pElemBlock: rgElemBlocks) {
                         const auto eElemType = pElemBlock->GetElementType();

                         const octave_idx_type iNumElem = pElemBlock->iGetNumElem();

                         if (!iNumElem) {
                              continue;
                         }

                         const octave_idx_type iNumColloc = pElemBlock->iGetNumCollocPoints(0, eMatType);

                         if (!iNumColloc) {
                              continue;
                         }

                         const dim_vector dimXg{iNumElem, iNumColloc, 3};

                         oSolution.SetField(PostProcData::VEC_EL_COLLOC_POINTS_RE, eElemType, NDArray{dimXg, 0.});

                         pElemBlock->PostProcElem(eMatType, oSolution, oParaOpt);

                         const NDArray& Xg = oSolution.GetField(PostProcData::VEC_EL_COLLOC_POINTS_RE, eElemType);
                         const char* pszElemName = ElementTypes::GetType(eElemType).name;

                         mapCollocPoints.assign(pszElemName, Xg);
                    }

                    retval.append(mapCollocPoints);
               } break;
               case Element::VEC_STRESS_CAUCH:
               case Element::VEC_STRAIN_TOTAL:
               case Element::SCA_STRESS_VMIS: {
                    if (nargin <= 4) {
                         throw std::runtime_error("fem_ass_matrix: argument sol is not optional for selected matrix type in argument matrix_type");
                    }

                    constexpr octave_idx_type iNumStress = 6;
                    const octave_idx_type iNumSteps = oSolution.GetNumSteps();

                    octave_scalar_map s_StressStrain, s_StressStrainm, s_vmis;

                    for (octave_idx_type j = 0; j < ElementTypes::iGetNumTypes(); ++j) {
                         const ElementTypes::TypeInfo& oElemType = ElementTypes::GetType(j);

                         if (!rgElemUse[oElemType.type]) {
                              continue;
                         }

                         switch (oElemType.type) {
                         case ElementTypes::ELEM_ISO8:
                         case ElementTypes::ELEM_ISO20:
                         case ElementTypes::ELEM_ISO20R:
                         case ElementTypes::ELEM_ISO27:
                         case ElementTypes::ELEM_PENTA15:
                         case ElementTypes::ELEM_TET10H:
                         case ElementTypes::ELEM_TET10:
                         case ElementTypes::ELEM_TET20: {
                              const auto iter_elem = elements.seek(oElemType.name);

                              if (iter_elem == elements.end()) {
                                   continue;
                              }

                              const int32NDArray elem_nodes = elements.contents(iter_elem).int32_array_value();

                              const octave_idx_type iNumElem = elem_nodes.rows();
                              const octave_idx_type iNumNodesElem = elem_nodes.columns();

                              const Element::FemMatrixType eMatTypeStressStrain = eMatType == Element::SCA_STRESS_VMIS ? Element::VEC_STRESS_CAUCH : eMatType;

                              PostProcData::FieldTypeReal eFieldType;

                              switch (eMatType) {
                              case Element::VEC_STRESS_CAUCH:
                              case Element::SCA_STRESS_VMIS:
                                   eFieldType = PostProcData::VEC_EL_STRUCT_STRESS_CAUCH_RE;
                                   break;
                              case Element::VEC_STRAIN_TOTAL:
                                   eFieldType = PostProcData::VEC_EL_STRUCT_STRAIN_TOTAL_RE;
                                   break;
                              default:
                                   FEM_ASSERT(0);
                                   throw std::logic_error("fem_ass_matrix: unexpected matrix type");
                              }

                              oSolution.SetField(eFieldType, oElemType.type, NDArray(dim_vector(iNumElem, iNumNodesElem, iNumStress, iNumSteps), 0.));

                              for (const auto& pElemBlock: rgElemBlocks) {
                                   if (pElemBlock->GetElementType() == oElemType.type) {
                                        pElemBlock->PostProcElem(eMatTypeStressStrain, oSolution, oParaOpt);
                                   }
                              }

                              const NDArray oElemStressStrain = oSolution.GetField(eFieldType, oElemType.type);

                              FEM_ASSERT(oElemStressStrain.ndims() == 4 || (oElemStressStrain.ndims() == 3 && iNumSteps == 1));
                              FEM_ASSERT(oElemStressStrain.dims()(0) == iNumElem);
                              FEM_ASSERT(oElemStressStrain.dims()(1) == iNumNodesElem);
                              FEM_ASSERT(oElemStressStrain.dims()(2) == iNumStress);
                              FEM_ASSERT(oElemStressStrain.ndims() > 3 ? oElemStressStrain.dims()(3) == iNumSteps : iNumSteps == 1);

                              Array<octave_idx_type> iStressStrain(dim_vector(iNumNodes, 1), 0);

                              for (octave_idx_type l = 0; l < iNumNodesElem; ++l) {
                                   for (octave_idx_type k = 0; k < iNumElem; ++k) {
                                        const octave_idx_type inode = elem_nodes.xelem(k + iNumElem * l).value() - 1;
                                        ++iStressStrain.xelem(inode);
                                   }
                              }

                              NDArray oNodalStressStrain(dim_vector(iNumNodes, iNumStress, iNumSteps), 0.);

                              for (octave_idx_type n = 0; n < iNumSteps; ++n) {
                                   for (octave_idx_type m = 0; m < iNumStress; ++m) {
                                        for (octave_idx_type l = 0; l < iNumNodesElem; ++l) {
                                             for (octave_idx_type k = 0; k < iNumElem; ++k) {
                                                  const octave_idx_type inode = elem_nodes.xelem(k + iNumElem * l).value() - 1;
                                                  oNodalStressStrain.xelem(inode + iNumNodes * (m + iNumStress *  n)) += oElemStressStrain.xelem(k + iNumElem * (l + iNumNodesElem * (m + n * iNumStress))) / iStressStrain.xelem(inode);
                                             }
                                        }
                                   }
                              }

                              NDArray oContStressStrain(oElemStressStrain.dims());

                              for (octave_idx_type n = 0; n < iNumSteps; ++n) {
                                   for (octave_idx_type m = 0; m < iNumStress; ++m) {
                                        for (octave_idx_type l = 0; l < iNumNodesElem; ++l) {
                                             for (octave_idx_type k = 0; k < iNumElem; ++k) {
                                                  const octave_idx_type inode = elem_nodes.xelem(k + iNumElem * l).value() - 1;
                                                  oContStressStrain.xelem(k + iNumElem * (l + iNumNodesElem * (m + n * iNumStress))) = oNodalStressStrain.xelem(inode + iNumNodes * (m + iNumStress * n));
                                             }
                                        }
                                   }
                              }

                              if (eMatType == Element::SCA_STRESS_VMIS) {
                                   NDArray vmis(dim_vector(iNumElem, iNumNodesElem, iNumSteps));

                                   for (octave_idx_type n = 0; n < iNumSteps; ++n) {
                                        for (octave_idx_type l = 0; l < iNumNodesElem; ++l) {
                                             for (octave_idx_type k = 0; k < iNumElem; ++k) {
                                                  const octave_idx_type ioffset = n * iNumStress;

                                                  const double tauxx = oContStressStrain.xelem(k + iNumElem * (l + iNumNodesElem * ioffset));
                                                  const double tauyy = oContStressStrain.xelem(k + iNumElem * (l + iNumNodesElem * (ioffset + 1)));
                                                  const double tauzz = oContStressStrain.xelem(k + iNumElem * (l + iNumNodesElem * (ioffset + 2)));
                                                  const double tauxy = oContStressStrain.xelem(k + iNumElem * (l + iNumNodesElem * (ioffset + 3)));
                                                  const double tauyz = oContStressStrain.xelem(k + iNumElem * (l + iNumNodesElem * (ioffset + 4)));
                                                  const double tauzx = oContStressStrain.xelem(k + iNumElem * (l + iNumNodesElem * (ioffset + 5)));

                                                  vmis.xelem(k + iNumElem * (l + iNumNodesElem * n)) = sqrt(tauxx * tauxx + tauyy * tauyy + tauzz * tauzz
                                                                                                            - (tauxx * tauyy + tauyy * tauzz + tauxx * tauzz)
                                                                                                            + 3. * (tauxy * tauxy + tauyz * tauyz + tauzx * tauzx));

                                             }
                                        }
                                   }

                                   s_vmis.assign(oElemType.name, vmis);
                              } else {
                                   s_StressStrain.assign(oElemType.name, oElemStressStrain);
                                   s_StressStrainm.assign(oElemType.name, oContStressStrain);
                              }
                         } break;

                         default:
                              break;
                         }
                    }

                    octave_scalar_map mapStressStrain;

                    switch (eMatType) {
                    case Element::SCA_STRESS_VMIS:
                         mapStressStrain.assign("vmis", s_vmis);
                         break;
                    case Element::VEC_STRESS_CAUCH:
                         mapStressStrain.assign("tau", s_StressStrain);
                         mapStressStrain.assign("taum", s_StressStrainm);
                         break;
                    case Element::VEC_STRAIN_TOTAL:
                         mapStressStrain.assign("epsilon", s_StressStrain);
                         mapStressStrain.assign("epsilonm", s_StressStrainm);
                         break;
                    default:
                         FEM_ASSERT(false);
                    }

                    retval.append(mapStressStrain);
               } break;
               case Element::VEC_PARTICLE_VELOCITY:
               case Element::VEC_PARTICLE_VELOCITY_C:
               case Element::SCA_ACOUSTIC_INTENSITY:
               case Element::SCA_ACOUSTIC_INTENSITY_C:
               {
                    if (nargin <= 4) {
                         throw std::runtime_error("fem_ass_matrix: argument sol is not optional for selected matrix type in argument matrix_type");
                    }

                    switch (eMatType) {
                    case Element::VEC_PARTICLE_VELOCITY:
                    case Element::SCA_ACOUSTIC_INTENSITY:
                         retval.append(AcousticPostProc<double>(rgElemUse, rgElemBlocks, elements, nodes, oSolution, eMatType, oParaOpt));
                         break;
                    case Element::VEC_PARTICLE_VELOCITY_C:
                    case Element::SCA_ACOUSTIC_INTENSITY_C:
                         retval.append(AcousticPostProc<std::complex<double> >(rgElemUse, rgElemBlocks, elements, nodes, oSolution, eMatType, oParaOpt));
                         break;
                    default:
                         FEM_ASSERT(false);
                    }
               } break;
               case Element::VEC_SURFACE_NORMAL_VECTOR:
                    retval.append(SurfaceNormalVectorPostProc(rgElemUse, rgElemBlocks, elements, nodes, oSolution, oParaOpt));
                    break;
               case Element::VEC_SURFACE_AREA:
                    retval.append(SurfaceAreaPostProc(rgElemUse, rgElemBlocks, nodes, oSolution, oParaOpt));
                    break;
               default:
                    throw std::runtime_error("fem_ass_matrix: invalid value for argument matrix_type");
               }
          }

          octave_scalar_map mat_info;
          ColumnVector beta(matrix_type.numel());

          for (octave_idx_type i = 0; i < matrix_type.numel(); ++i) {
               auto eMatType = static_cast<Element::FemMatrixType>(matrix_type.xelem(i).value());
               beta.xelem(i) = oMatAss.bHaveMatrixInfo(eMatType) ? oMatAss.GetMatrixInfo(eMatType).beta : 0.;
          }

          mat_info.assign("beta", beta);

          retval.append(mat_info);
          retval.append(oMeshInfo.Get());

          double detJmin = oMeshInfo.dGet(MeshInfo::JACOBIAN_DET, MeshInfo::STAT_MIN);

          if (detJmin <= 0.) {
               warning_with_id("mboct-fem-pkg:invalid-mesh", "Jacobian of solid element is singular or negative det(J)=%g", detJmin);
          }
     } catch (const std::exception& err) {
          error("%s", err.what());
          return retval;
     }

     return retval;
}

#define DEFINE_GLOBAL_CONSTANT(NAMESPACE, CONST, DESCRIPTION)   \
     DEFUN_DLD(FEM_##CONST, args, nargout,                      \
               "-*- texinfo -*-\n"                              \
               "@deftypefn {} @var{id} = FEM_" #CONST  "()\n"   \
               DESCRIPTION "\n"                                 \
               "@end deftypefn\n")                              \
                                                                \
     {                                                          \
     return octave_value(octave_int32(NAMESPACE::CONST));       \
     }

DEFINE_GLOBAL_CONSTANT(Element, MAT_DAMPING, "complete damping matrix")
DEFINE_GLOBAL_CONSTANT(Element, MAT_DAMPING_SYM, "upper triangular part of the damping matrix")
DEFINE_GLOBAL_CONSTANT(Element, MAT_DAMPING_SYM_L, "lower triangular part of the damping matrix")
DEFINE_GLOBAL_CONSTANT(Element, MAT_INERTIA_INV3, "MBDyn's invariant 3 with respect to the global reference frame")
DEFINE_GLOBAL_CONSTANT(Element, MAT_INERTIA_INV4, "MBDyn's invariant 4 with respect to the global reference frame")
DEFINE_GLOBAL_CONSTANT(Element, MAT_INERTIA_INV5, "MBDyn's invarinat 5 with respect to the global reference frame")
DEFINE_GLOBAL_CONSTANT(Element, MAT_INERTIA_INV8, "MBDyn's invariant 8 with respect to the global reference frame")
DEFINE_GLOBAL_CONSTANT(Element, MAT_INERTIA_INV9, "MBDyn's invariant 9 with respect to the global reference frame")
DEFINE_GLOBAL_CONSTANT(Element, MAT_INERTIA_J, "3x3 inertia matrix with respect to the origin of the global reference frame")
DEFINE_GLOBAL_CONSTANT(Element, MAT_MASS, "complete consistent mass matrix")
DEFINE_GLOBAL_CONSTANT(Element, MAT_MASS_LUMPED, "diagonal lumped mass matrix")
DEFINE_GLOBAL_CONSTANT(Element, MAT_MASS_SYM, "upper triangular part of the mass matrix")
DEFINE_GLOBAL_CONSTANT(Element, MAT_MASS_SYM_L, "lower triangular part of the mass matrix")
DEFINE_GLOBAL_CONSTANT(Element, MAT_STIFFNESS, "complete stiffness matrix")
DEFINE_GLOBAL_CONSTANT(Element, MAT_STIFFNESS_IM, "imaginary part of complete stiffness matrix")
DEFINE_GLOBAL_CONSTANT(Element, MAT_STIFFNESS_SYM, "upper triangular part of the stiffness matrix")
DEFINE_GLOBAL_CONSTANT(Element, MAT_STIFFNESS_SYM_L, "lower triangular part of the stiffness matrix")
DEFINE_GLOBAL_CONSTANT(Element, MAT_STIFFNESS_TAU0, "stiffness matrix resulting from pre-stress")
DEFINE_GLOBAL_CONSTANT(Element, MAT_STIFFNESS_OMEGA, "stiffness matrix resulting from a rotating reference frame with constant angular velocity")
DEFINE_GLOBAL_CONSTANT(Element, MAT_STIFFNESS_OMEGA_DOT, "stiffness matrix resulting from a rotating reference frame with constant angular acceleration")
DEFINE_GLOBAL_CONSTANT(Element, MAT_DAMPING_OMEGA, "damping matrix resulting from a rotating reference frame")
DEFINE_GLOBAL_CONSTANT(Element, SCA_STRESS_VMIS, "Van Mises stress")
DEFINE_GLOBAL_CONSTANT(Element, SCA_TOT_MASS, "total mass of all elements")
DEFINE_GLOBAL_CONSTANT(Element, VEC_INERTIA_M1, "center of mass times total mass of all elements with respect to the origin of the global reference frame")
DEFINE_GLOBAL_CONSTANT(Element, VEC_LOAD_CONSISTENT, "consistent load vector")
DEFINE_GLOBAL_CONSTANT(Element, VEC_LOAD_LUMPED, "lumped load vector")
DEFINE_GLOBAL_CONSTANT(Element, VEC_STRESS_CAUCH, "linear stress tensor")
DEFINE_GLOBAL_CONSTANT(Element, VEC_STRAIN_TOTAL, "linear strain tensor")
DEFINE_GLOBAL_CONSTANT(SurfToNodeConstrBase, CT_FIXED, "build constraints in all three directions in space")
DEFINE_GLOBAL_CONSTANT(SurfToNodeConstrBase, CT_SLIDING, "build only one constraint normal to the surface")
DEFINE_GLOBAL_CONSTANT(Element, MAT_THERMAL_COND, "thermal conductivity matrix")
DEFINE_GLOBAL_CONSTANT(Element, MAT_HEAT_CAPACITY, "heat capacity matrix")
DEFINE_GLOBAL_CONSTANT(Element, VEC_LOAD_THERMAL, "thermal load vector")
DEFINE_GLOBAL_CONSTANT(Element, MAT_STIFFNESS_ACOUSTICS_RE, "real part of acoustic stiffness matrix")
DEFINE_GLOBAL_CONSTANT(Element, MAT_STIFFNESS_ACOUSTICS_IM, "imaginary part of acoustic stiffness matrix")
DEFINE_GLOBAL_CONSTANT(Element, MAT_MASS_ACOUSTICS_RE, "real acoustic mass matrix")
DEFINE_GLOBAL_CONSTANT(Element, MAT_MASS_ACOUSTICS_IM, "imaginary acoustic mass matrix")
DEFINE_GLOBAL_CONSTANT(Element, MAT_DAMPING_ACOUSTICS_RE, "real part of acoustic damping matrix")
DEFINE_GLOBAL_CONSTANT(Element, MAT_DAMPING_ACOUSTICS_IM, "imaginary part of acoustic damping matrix")
DEFINE_GLOBAL_CONSTANT(Element, VEC_LOAD_ACOUSTICS, "acoustic load vector")
DEFINE_GLOBAL_CONSTANT(Element, VEC_PARTICLE_VELOCITY, "acoustic particle velocity")
DEFINE_GLOBAL_CONSTANT(Element, SCA_ACOUSTIC_INTENSITY, "acoustic intensity and sound power")
DEFINE_GLOBAL_CONSTANT(Element, VEC_PARTICLE_VELOCITY_C, "complex acoustic particle velocity")
DEFINE_GLOBAL_CONSTANT(Element, SCA_ACOUSTIC_INTENSITY_C, "acoustic intensity and sound power for complex solutions")
DEFINE_GLOBAL_CONSTANT(Element, VEC_SURFACE_NORMAL_VECTOR, "surface normal vector at elements")
DEFINE_GLOBAL_CONSTANT(Element, VEC_SURFACE_AREA, "surface area for pressure loads")
DEFINE_GLOBAL_CONSTANT(Element, MAT_STIFFNESS_FLUID_STRUCT_RE, "real fluid-structure interaction stiffness matrix")
DEFINE_GLOBAL_CONSTANT(Element, MAT_STIFFNESS_FLUID_STRUCT_IM, "imaginary fluid-structure interaction stiffness matrix")
DEFINE_GLOBAL_CONSTANT(Element, MAT_MASS_FLUID_STRUCT_RE, "real fluid-structure interaction mass matrix")
DEFINE_GLOBAL_CONSTANT(Element, MAT_MASS_FLUID_STRUCT_IM, "imaginary fluid-structure interaction mass matrix")
DEFINE_GLOBAL_CONSTANT(Element, MAT_DAMPING_FLUID_STRUCT_RE, "real part of fluid-structure interaction damping matrix")
DEFINE_GLOBAL_CONSTANT(Element, MAT_DAMPING_FLUID_STRUCT_IM, "imaginary part of fluid-structure interaction damping matrix")
DEFINE_GLOBAL_CONSTANT(Element, VEC_LOAD_FLUID_STRUCT, "fluid-structure interaction load vector")
DEFINE_GLOBAL_CONSTANT(Element, VEC_COLL_MASS, "collocation points of mass matrix")
DEFINE_GLOBAL_CONSTANT(Element, VEC_COLL_STIFFNESS, "collocation points of stiffness matrix")
DEFINE_GLOBAL_CONSTANT(Element, VEC_COLL_HEAT_CAPACITY, "collocation points of heat capacity matrix")
DEFINE_GLOBAL_CONSTANT(Element, VEC_COLL_THERMAL_COND, "collocation points of thermal conductivity matrix")
DEFINE_GLOBAL_CONSTANT(Element, VEC_COLL_MASS_ACOUSTICS, "collocation points of acoustic mass matrix")
DEFINE_GLOBAL_CONSTANT(Element, VEC_COLL_STIFF_ACOUSTICS, "collocation points of acoustic stiffness matrix")
DEFINE_GLOBAL_CONSTANT(Element, VEC_COLL_MASS_FLUID_STRUCT, "collocation points of fluid structure interaction mass matrix")
DEFINE_GLOBAL_CONSTANT(Element, VEC_COLL_STIFF_FLUID_STRUCT, "collocation points of fluid structure interaction stiffness matrix")
DEFINE_GLOBAL_CONSTANT(DofMap, DO_STRUCTURAL, "structural domain")
DEFINE_GLOBAL_CONSTANT(DofMap, DO_THERMAL, "thermal domain")
DEFINE_GLOBAL_CONSTANT(DofMap, DO_ACOUSTICS, "acoustic domain")
DEFINE_GLOBAL_CONSTANT(DofMap, DO_FLUID_STRUCT, "fluid-structure interaction domain")
