// Copyright (C) 2018(-2021) Reinhard <octave-user@a1.net>
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
#include <array>
#include <cassert>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
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

class DofMap {
public:
     enum DomainType: unsigned {
          DO_STRUCTURAL = 0x00000000u,
          DO_THERMAL    = 0x01000000u,
     };

     enum ElementType {
          ELEM_RBE3 = 0,
          ELEM_JOINT,
          ELEM_TYPE_COUNT,
          ELEM_NODOF = -1
     };

     explicit DofMap(DomainType eDomain, const int32NDArray& ndof, const array<int32NDArray, ELEM_TYPE_COUNT>& edof, octave_idx_type totdof)
          :eDomain(eDomain), ndof(ndof), edof(edof), totdof(totdof) {
     }

     explicit DofMap(DomainType eDomain, const int32NDArray& ndof, octave_idx_type totdof)
          :eDomain(eDomain), ndof(ndof), totdof(totdof) {
     }

     octave_idx_type GetNodeDofIndex(octave_idx_type inode, octave_idx_type idof) const {
          return ndof.xelem(inode, idof);
     }

     octave_idx_type GetElemDofIndex(ElementType eElemType, octave_idx_type ielem, octave_idx_type idof) const {
          return edof[eElemType].xelem(ielem, idof).value();
     }

     octave_idx_type iGetNumDof() const {
          return totdof;
     }

     octave_idx_type iGetNodeMaxDofIndex() const {
          return ndof.columns();
     }

     static octave_idx_type iGetNodeMaxDofIndex(DomainType eDomain) {
          switch (eDomain) {
          case DO_STRUCTURAL:
               return 6;
               
          case DO_THERMAL:
               return 1;
               
          default:
               return -1;
          }
     }

     DomainType GetDomain() const {
          return eDomain;
     }

private:
     DomainType eDomain;
     int32NDArray ndof;
     array<int32NDArray, ELEM_TYPE_COUNT> edof;
     octave_idx_type totdof;
};

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

          static const char szStatName[STAT_COUNT][5] = {
               "min",
               "max",
               "mean"
          };

          static const char szInfoName[INFO_COUNT][6] = {
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
     Material(const Matrix& C, double rho, double alpha, double beta, double gamma, const Matrix& k, double cp)
          :rho(rho),
           alpha(alpha),
           beta(beta),
           gamma(gamma),
           cp(cp),
           C(C), k(k) {
          FEM_ASSERT(C.rows() == 6);
          FEM_ASSERT(C.columns() == 6);

          const double c = C.xelem(0, 0);
          const double b = C.xelem(5, 5) / c;

          E = ((4 * b * b - 3 * b) * c) / (b - 1);
          nu = (2 * b - 1) / (2 * b - 2);
     }

     Material(double E, double nu, double rho, double alpha, double beta, double gamma, const Matrix& k, double cp)
          :E(E), nu(nu), rho(rho), alpha(alpha), beta(beta), gamma(gamma), cp(cp), C(6, 6, 0.), k(k) {

          const double a = nu / (1 - nu);
          const double b = (1 - 2 * nu) / (2 * (1 - nu));
          const double c = E * (1 - nu) / ((1 + nu) * (1 - 2 * nu));

          for (octave_idx_type j = 0; j < 3; ++j) {
               for (octave_idx_type i = 0; i < 3; ++i) {
                    C.xelem(i, j) = (i == j) ? c : a * c;
               }

               C.xelem(j + 3, j + 3) = b * c;
          }
     }

     Material(const Material& oMat)
          :rho(oMat.rho), alpha(oMat.alpha), beta(oMat.beta), gamma(oMat.gamma), cp(oMat.cp), C(oMat.C), k(oMat.k) {
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

private:
     double E, nu, rho, alpha, beta, gamma, cp;
     Matrix C, k;
};

class IntegrationRule
{
public:
     IntegrationRule() {
     }

     void SetNumEvalPoints(octave_idx_type iNumEvalPoints, octave_idx_type iNumDirections) {
          r.resize(iNumEvalPoints, iNumDirections);
          alpha.resize(iNumEvalPoints);
     }

     void SetPosition(octave_idx_type iEvalPnt, octave_idx_type iDirection, double ri) {
          r.xelem(iEvalPnt, iDirection) = ri;
     }

     void SetWeight(octave_idx_type iEvalPnt, double alphai) {
          alpha.xelem(iEvalPnt) = alphai;
     }

     double dGetPosition(octave_idx_type iEvalPnt, octave_idx_type iDirection) const {
          return r.xelem(iEvalPnt, iDirection);
     }

     double dGetWeight(octave_idx_type iEvalPnt) const {
          return alpha.xelem(iEvalPnt);
     }

     octave_idx_type iGetNumDirections() const {
          return r.columns();
     }

     octave_idx_type iGetNumEvalPoints() const {
          FEM_ASSERT(r.rows() == alpha.rows());

          return r.rows();
     }

private:
     Matrix r;
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
                    throw std::runtime_error("argument load_case.dTheta must be a real column vector with the same number of rows like mesh.nodes");
               }
          }

          for (octave_idx_type i = 0; i < rgRefStrain.numel(); ++i) {
               const octave_value ov_RefStrain = rgRefStrain.xelem(i);

               if (!(ov_RefStrain.isstruct() && ov_RefStrain.numel() == 1)) {
                    throw std::runtime_error("argument load_case.epsilon0 must be a scalar struct");
               }
          }          
     }

     Cell rgTemperature;
     Cell rgRefStrain;
};

class ElementTypes {
public:
     enum TypeId {
          ELEM_ISO8 = 0,
          ELEM_ISO20,
          ELEM_PENTA15,
          ELEM_TET10H,
          ELEM_TET10,
          ELEM_BEAM2,
          ELEM_RBE3,
          ELEM_JOINT,
          ELEM_SFNCON4,
          ELEM_SFNCON6,
          ELEM_SFNCON6H,
          ELEM_SFNCON8,
          ELEM_PRESSURE_ISO4,
          ELEM_PRESSURE_QUAD8,
          ELEM_PRESSURE_TRIA6,
          ELEM_PRESSURE_TRIA6H,
          ELEM_STRUCT_FORCE,
          ELEM_THERM_CONV_ISO4,
          ELEM_THERM_CONV_QUAD8,
          ELEM_THERM_CONV_TRIA6,
          ELEM_THERM_CONV_TRIA6H,
          ELEM_TYPE_COUNT
     };

     struct TypeInfo {
          TypeId type;
          char name[9];
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
     static const TypeInfo rgElemTypes[ELEM_TYPE_COUNT];
};

const ElementTypes::TypeInfo ElementTypes::rgElemTypes[ElementTypes::ELEM_TYPE_COUNT] = {
     {ElementTypes::ELEM_ISO8,             "iso8",     8,  8, DofMap::ELEM_NODOF},
     {ElementTypes::ELEM_ISO20,            "iso20",   20, 20, DofMap::ELEM_NODOF},
     {ElementTypes::ELEM_PENTA15,          "penta15", 15, 15, DofMap::ELEM_NODOF},
     {ElementTypes::ELEM_TET10H,           "tet10h",  10, 10, DofMap::ELEM_NODOF},
     {ElementTypes::ELEM_TET10,            "tet10",   10, 10, DofMap::ELEM_NODOF},
     {ElementTypes::ELEM_BEAM2,            "beam2",    2,  2, DofMap::ELEM_NODOF},
     {ElementTypes::ELEM_RBE3,             "rbe3",     2, -1, DofMap::ELEM_RBE3},
     {ElementTypes::ELEM_JOINT,            "joints",   1, -1, DofMap::ELEM_JOINT},
     {ElementTypes::ELEM_SFNCON4,          "sfncon4",  1, -1, DofMap::ELEM_JOINT},
     {ElementTypes::ELEM_SFNCON6,          "sfncon6",  1, -1, DofMap::ELEM_JOINT},
     {ElementTypes::ELEM_SFNCON6H,         "sfncon6h", 1, -1, DofMap::ELEM_JOINT},
     {ElementTypes::ELEM_SFNCON8,          "sfncon8",  1, -1, DofMap::ELEM_JOINT},
     {ElementTypes::ELEM_PRESSURE_ISO4,    "iso4",     4,  4, DofMap::ELEM_NODOF},
     {ElementTypes::ELEM_PRESSURE_QUAD8,   "quad8",    8,  8, DofMap::ELEM_NODOF},
     {ElementTypes::ELEM_PRESSURE_TRIA6,   "tria6",    6,  6, DofMap::ELEM_NODOF},
     {ElementTypes::ELEM_PRESSURE_TRIA6H,  "tria6h",   6,  6, DofMap::ELEM_NODOF},     
     {ElementTypes::ELEM_STRUCT_FORCE,     "force",    1, -1, DofMap::ELEM_NODOF},
     {ElementTypes::ELEM_THERM_CONV_ISO4,  "iso4",     4,  4, DofMap::ELEM_NODOF},
     {ElementTypes::ELEM_THERM_CONV_QUAD8, "quad8",    8,  8, DofMap::ELEM_NODOF},
     {ElementTypes::ELEM_THERM_CONV_TRIA6, "tria6",    6,  6, DofMap::ELEM_NODOF},
     {ElementTypes::ELEM_THERM_CONV_TRIA6H,"tria6h",   6,  6, DofMap::ELEM_NODOF}
};

class Element
{
public:
     enum MatSymmetryFlag: unsigned {
          MAT_SYM_FULL = 0x0u,
          MAT_SYM_UPPER = 0x1u,
          MAT_SYM_LOWER = 0x2u,
          MAT_SYM_DIAG = 0x3u,
          MAT_SYM_MASK = 0xFu
     };

     enum FemMatrixType: unsigned {
          MAT_UNKNOWN = static_cast<unsigned>(-1),
          MAT_STIFFNESS = 0x100u,
          MAT_STIFFNESS_SYM = MAT_STIFFNESS | MAT_SYM_UPPER,
          MAT_STIFFNESS_SYM_L = MAT_STIFFNESS | MAT_SYM_LOWER,
          MAT_MASS = 0x200u,
          MAT_MASS_SYM = MAT_MASS | MAT_SYM_UPPER,
          MAT_MASS_SYM_L = MAT_MASS | MAT_SYM_LOWER,
          MAT_MASS_LUMPED = MAT_MASS | MAT_SYM_DIAG,
          MAT_DAMPING = 0x300u,
          MAT_DAMPING_SYM = MAT_DAMPING | MAT_SYM_UPPER,
          MAT_DAMPING_SYM_L = MAT_DAMPING | MAT_SYM_LOWER,
          SCA_TOT_MASS = 0x400u,
          VEC_INERTIA_M1 = 0x500u,
          MAT_INERTIA_J = 0x600u,
          MAT_INERTIA_INV3 = 0x700u,
          MAT_INERTIA_INV4 = 0x800u,
          MAT_INERTIA_INV5 = 0x900u,
          MAT_INERTIA_INV8 = 0xA00u,
          MAT_INERTIA_INV9 = 0xB00u,
          MAT_ACCEL_LOAD = 0xC00u,
          VEC_LOAD_CONSISTENT = 0xD00u,
          VEC_LOAD_LUMPED = 0xE00u,
          VEC_STRESS_CAUCH = 0xF00u,
          VEC_STRAIN_TOTAL = 0x1000u,
          SCA_STRESS_VMIS = 0x2000u,
          MAT_THERMAL_COND = 0x4000u | DofMap::DO_THERMAL,
          MAT_HEAT_CAPACITY = 0x8000u | DofMap::DO_THERMAL,
          VEC_LOAD_THERMAL = 0x10000u | DofMap::DO_THERMAL
     };

     Element(ElementTypes::TypeId eltype, octave_idx_type id, const Matrix& X, const Material* material, const int32NDArray& nodes)
          :eltype(eltype), id(id), X(X), material(material), nodes(nodes) {

          FEM_ASSERT(X.columns() == nodes.numel());
     }

     Element(const Element& oElem)
          :eltype(oElem.eltype), id(oElem.id), X(oElem.X), material(oElem.material), nodes(oElem.nodes) {
     }

     virtual ~Element() {
     }

     virtual void Assemble(MatrixAss& mat, MeshInfo& info, const DofMap& dof, FemMatrixType eMatType) const=0;

     virtual void PostProcElem(NDArray& mat, FemMatrixType eMatType, const NDArray& U) const {
     }

     virtual octave_idx_type iGetWorkSpaceSize(FemMatrixType eMatType) const=0;

     virtual double dGetMass() const {
          return 0.;
     }

     virtual bool bNeedMatrixInfo(FemMatrixType eMatType) const {
          return false;
     }

     virtual void Extract(octave_idx_type& idx, octave_map& sElem) const {
     }

protected:
     const ElementTypes::TypeId eltype;
     const octave_idx_type id;
     const Matrix X;
     const Material* material;
     const int32NDArray nodes;
};

class MatrixAss {
public:
     struct MatrixInfo {
          MatrixInfo()
               :alpha(1.), beta(1.) {
          }

          double alpha; // scale factor for double Lagrange multipliers
          double beta; // scale factor for constraint equations
     };

     explicit MatrixAss(octave_idx_type max_nnz)
          :eMatType(Element::MAT_UNKNOWN),
           nnz(0),
           ridx(dim_vector(max_nnz, 1), 0),
           cidx(dim_vector(max_nnz, 1), 0),
           data(max_nnz, 0.) {
     }

     const MatrixInfo& GetMatrixInfo() const {
          return info;
     }

     void UpdateMatrixInfo() {
          ColumnVector diagA(nnz, 0.);

          for (octave_idx_type i = 0; i < nnz; ++i) {
               if (ridx.xelem(i).value() == cidx.xelem(i).value()) {
                    diagA.xelem(ridx.xelem(i)) += data.xelem(i);
               }
          }

          constexpr double INIT_MIN = std::numeric_limits<double>::max();
          constexpr double INIT_MAX = -std::numeric_limits<double>::max();

          double minA = INIT_MIN;
          double maxA = INIT_MAX;

          for (octave_idx_type i = 0; i < nnz; ++i) {
               const double absA = fabs(diagA.xelem(i));

               if (absA > maxA) {
                    maxA = absA;
               }

               if (absA < minA) {
                    minA = absA;
               }
          }

          // According to Code_Aster r3.03.08
          info.alpha = info.beta = (minA != INIT_MIN && maxA != INIT_MAX) ? 0.5 * (minA + maxA) : 1.;
     }

     void Insert(double d, octave_idx_type r, octave_idx_type c) {
          if (bNeedToInsertElem(r, c)) {
               Resize(nnz + 1);
               InsertRaw(d, r, c);
          }
     }

     void Insert(const Matrix& Ke, const int32NDArray& r, const int32NDArray& c) {
          for (octave_idx_type j = 0; j < Ke.columns(); ++j) {
               for (octave_idx_type i = 0; i < Ke.rows(); ++i) {
                    Insert(Ke.xelem(i, j), r.xelem(i), c.xelem(j));
               }
          }
     }

     void Finish() {
          // Do not resize the workspace here because it could be reused for other matrices!
     }

     const int32NDArray& RowIndex() const { return ridx; }
     const int32NDArray& ColIndex() const { return cidx; }
     const ColumnVector& Data() const { return data; }

     SparseMatrix Assemble(const DofMap& oDofMap, octave_idx_type iNumLoads) const {
          FEM_ASSERT(eMatType != Element::MAT_UNKNOWN);

          const octave_idx_type iNumRows = oDofMap.iGetNumDof();
          octave_idx_type iNumCols = -1;

          switch (eMatType) {
          case Element::VEC_LOAD_CONSISTENT:
          case Element::VEC_LOAD_LUMPED:
          case Element::VEC_LOAD_THERMAL:
               iNumCols = iNumLoads;
               break;

          case Element::MAT_ACCEL_LOAD:
               iNumCols = 3;
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

     void Reset(Element::FemMatrixType eMatTypeCurr, const MatrixInfo& oMatInfo) {
          eMatType = eMatTypeCurr;
          nnz = 0;
          info = oMatInfo;
     }

private:
     void Resize(octave_idx_type nnz_new) {
          if (data.rows() < nnz_new) {
#ifdef DEBUG
               throw std::runtime_error("fem_ass_matrix: allocated workspace size exceeded");
#endif
               data.resize(nnz_new, 0.);
               ridx.resize(dim_vector(nnz_new, 1), 0);
               cidx.resize(dim_vector(nnz_new, 1), 0);
          }
     }


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

          data.xelem(nnz) = d;
          ridx.xelem(nnz) = r;
          cidx.xelem(nnz) = c;
          ++nnz;
     }

     Element::FemMatrixType eMatType;
     octave_idx_type nnz;
     int32NDArray ridx, cidx;
     ColumnVector data;
     MatrixInfo info;
};

class ElemJoint: public Element
{
public:
     ElemJoint(ElementTypes::TypeId eltype, octave_idx_type id, const Matrix& X, const Material* material, const int32NDArray& nodes, const Matrix& C, const Matrix& U)
          :Element(eltype, id, X, material, nodes), C(C), U(U) {
          FEM_ASSERT(C.columns() == nodes.numel() * 6);
          FEM_ASSERT(C.rows() <= C.columns());
          FEM_ASSERT(C.rows() >= 1);
          FEM_ASSERT(U.rows() == C.rows());
     }

     ElemJoint(const ElemJoint& oElem)
          :Element(oElem), C(oElem.C), U(oElem.U) {
     }

     virtual void Assemble(MatrixAss& mat, MeshInfo& info, const DofMap& dof, const FemMatrixType eMatType) const {
          switch (eMatType) {
          case MAT_STIFFNESS:
          case MAT_STIFFNESS_SYM:
          case MAT_STIFFNESS_SYM_L: {
               int32NDArray ndofidx(dim_vector(nodes.numel() * 6, 1), -1);

               for (octave_idx_type inode = 0; inode < nodes.numel(); ++inode) {
                    for (octave_idx_type idof = 0; idof < 6; ++idof) {
                         ndofidx.xelem(inode * 6 + idof) = dof.GetNodeDofIndex(nodes.xelem(inode).value() - 1, idof);
                    }
               }

               int32NDArray edofidx(dim_vector(C.rows(), 1));

               for (octave_idx_type idof = 0; idof < edofidx.numel(); ++idof) {
                    edofidx.xelem(idof) = dof.GetElemDofIndex(DofMap::ELEM_JOINT, id - 1, idof);
               }

               const double beta = mat.GetMatrixInfo().beta;

               for (octave_idx_type j = 0; j < C.columns(); ++j) {
                    for (octave_idx_type i = 0; i < C.rows(); ++i) {
                         const double Cij = beta * C.xelem(i, j);
                         mat.Insert(Cij, ndofidx.xelem(j), edofidx.xelem(i));
                         mat.Insert(Cij, edofidx.xelem(i), ndofidx.xelem(j));
                    }
               }
          } break;
          case VEC_LOAD_CONSISTENT:
          case VEC_LOAD_LUMPED: {
               int32NDArray edofidx(dim_vector(C.rows(), 1));

               for (octave_idx_type idof = 0; idof < edofidx.numel(); ++idof) {
                    edofidx.xelem(idof) = dof.GetElemDofIndex(DofMap::ELEM_JOINT, id - 1, idof);
               }

               const double beta = mat.GetMatrixInfo().beta;

               for (octave_idx_type j = 0; j < U.columns(); ++j) {
                    for (octave_idx_type i = 0; i < U.rows(); ++i) {
                         mat.Insert(beta * U.xelem(i, j), edofidx.xelem(i), j + 1);
                    }
               }
          } break;
          default:
               break;
          }
     }

     virtual octave_idx_type iGetWorkSpaceSize(FemMatrixType eMatType) const {
          switch (eMatType) {
          case MAT_STIFFNESS:
               return 2 * C.rows() * C.columns();
          case MAT_STIFFNESS_SYM:
          case MAT_STIFFNESS_SYM_L:
               return C.rows() * C.columns();
          case VEC_LOAD_CONSISTENT:
          case VEC_LOAD_LUMPED:
               return U.rows() * U.columns();
          default:
               return 0;
          }
     }

     virtual bool bNeedMatrixInfo(Element::FemMatrixType eMatType) const {
          switch (eMatType) {
          case MAT_STIFFNESS:
          case MAT_STIFFNESS_SYM:
          case MAT_STIFFNESS_SYM_L:
               return true;
          default:
               return false;
          }
     }

     virtual void Extract(octave_idx_type& idx, octave_map& sElem) const {
          FEM_ASSERT(sElem.numel() > idx);

          Cell& ovC = sElem.contents("C");
          Cell& ovNodes = sElem.contents("nodes");

          ovC(idx) = C;
          ovNodes(idx) = nodes.transpose();
          ++idx;
     }
private:
     const Matrix C, U;
};

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

     virtual void Assemble(MatrixAss& mat, MeshInfo& info, const DofMap& dof, const FemMatrixType eMatType) const {
          switch (eMatType) {
          case MAT_STIFFNESS:
          case MAT_STIFFNESS_SYM:
          case MAT_STIFFNESS_SYM_L:
               break;

          default:
               return;
          }

          int32NDArray ndofidx(dim_vector(nodes.numel() * 6, 1), -1);

          for (octave_idx_type inode = 0; inode < nodes.numel(); ++inode) {
               for (octave_idx_type idof = 0; idof < 6; ++idof) {
                    ndofidx.xelem(inode * 6 + idof) = dof.GetNodeDofIndex(nodes.xelem(inode).value() - 1, idof);
               }
          }

          int32NDArray edofidx(dim_vector(6, 1));

          for (octave_idx_type idof = 0; idof < edofidx.rows(); ++idof) {
               edofidx.xelem(idof) = dof.GetElemDofIndex(DofMap::ELEM_RBE3, id - 1, idof);
          }

          Matrix xi(3, nodes.numel() - 1);

          for (octave_idx_type j = 1; j < nodes.numel(); ++j) {
               for (octave_idx_type i = 0; i < 3; ++i) {
                    xi.xelem(i, j - 1) = X.xelem(i, j) - X.xelem(i, 0);
               }
          }

          FEM_TRACE("xi=[\n" << xi << "];\n");

          Matrix S(xi.columns() * 6, 6);

          for (octave_idx_type k = 0; k < xi.columns(); ++k) {
               for (octave_idx_type j = 0; j < 6; ++j) {
                    for (octave_idx_type i = 0; i < 6; ++i) {
                         const bool alpha = j < 3 || ndofidx.xelem((k + 1) * 6 + j).value() >= 0;
                         S.xelem(6 * k + i, j) = alpha * (i == j);
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

               S.xelem(6 * k + 0, 4) =  xi.xelem(2, k);
               S.xelem(6 * k + 0, 5) = -xi.xelem(1, k);
               S.xelem(6 * k + 1, 3) = -xi.xelem(2, k);
               S.xelem(6 * k + 1, 5) =  xi.xelem(0, k);
               S.xelem(6 * k + 2, 3) =  xi.xelem(1, k);
               S.xelem(6 * k + 2, 4) = -xi.xelem(0, k);
          }

          FEM_TRACE("S=[\n" << S << "];\n");

          double Lc2 = 0.;

          for (octave_idx_type k = 0; k < xi.columns(); ++k) {
               double norm_xik = 0.;

               for (octave_idx_type i = 0; i < 3; ++i) {
                    norm_xik += xi.xelem(i, k) * xi.xelem(i, k);
               }

               Lc2 += sqrt(norm_xik);
          }

          Lc2 /= xi.columns();
          Lc2 *= Lc2;

          FEM_TRACE("Lc2=" << Lc2 << ";\n");

          ColumnVector W(xi.columns() * 6);

          for (octave_idx_type k = 0; k < xi.columns(); ++k) {
               const double omegak = omega.xelem(k);

               for (octave_idx_type i = 0; i < 6; ++i) {
                    W.xelem(k * 6 + i) = omegak * (i < 3 ? 1. : Lc2);
               }
          }

          Matrix STWS(6, 6, 0.);

          for (octave_idx_type j = 0; j < 6; ++j) {
               for (octave_idx_type i = 0; i < 6; ++i) {
                    for (octave_idx_type k = 0; k < S.rows(); ++k) {
                         STWS.xelem(i, j) += S.xelem(k, i) * W.xelem(k) * S.xelem(k, j);
                    }
               }
          }

          FEM_TRACE("f(S^T*W*S)=[\n" << (S.transpose() * DiagMatrix(W) * S - STWS) << "];\n");

          octave_idx_type info;

          Matrix X(STWS.inverse(info));

          FEM_TRACE("X=[\n" << X << "];\n");

          if (info != 0) {
               std::ostringstream os;

               os << "rbe3 element id " << id << ": X matrix is singular";

               throw std::runtime_error(os.str());
          }

          Matrix B(xi.columns() * 6, 6);

          for (octave_idx_type l = 0; l < xi.columns(); ++l) {
               for (octave_idx_type j = 0; j < 6; ++j) {
                    for (octave_idx_type i = 0; i < 6; ++i) {
                         double Bijl = 0.;

                         for (octave_idx_type k = 0; k < 6; ++k) {
                              Bijl += W.xelem(l * 6 + i) * S.xelem(l * 6 + i, k) * X.xelem(k, j);
                         }

                         B.xelem(l * 6 + i, j) = Bijl;
                    }
               }
          }

          FEM_TRACE("f=[\n" << (B - DiagMatrix(W) * S * X) << "];\n");
          FEM_TRACE("W=[\n" << W << "];\n");
          FEM_TRACE("B=[\n" << B << "];\n");

          const double beta = mat.GetMatrixInfo().beta;

          for (octave_idx_type i = 0; i < 6; ++i) {
               mat.Insert(-beta, ndofidx.xelem(i), edofidx.xelem(i));
               mat.Insert(-beta, edofidx.xelem(i), ndofidx.xelem(i));
          }

          for (octave_idx_type j = 0; j < 6; ++j) {
               for (octave_idx_type i = 0; i < xi.columns() * 6; ++i) {
                    const double Bij = beta * B.xelem(i, j);
                    mat.Insert(Bij, ndofidx.xelem(i + 6), edofidx.xelem(j));
                    mat.Insert(Bij, edofidx.xelem(j), ndofidx.xelem(i + 6));
               }
          }
     }

     virtual octave_idx_type iGetWorkSpaceSize(FemMatrixType eMatType) const {
          switch (eMatType) {
          case MAT_STIFFNESS:
               return 8 * 6 + 4 * 6 * 6 * (X.columns() - 1);
          case MAT_STIFFNESS_SYM:
          case MAT_STIFFNESS_SYM_L:
               return 4 * 6 + 2 * 6 * 6 * (X.columns() - 1);
          default:
               return 0;
          }
     }

     virtual bool bNeedMatrixInfo(Element::FemMatrixType eMatType) const {
          switch (eMatType) {
          case MAT_STIFFNESS:
          case MAT_STIFFNESS_SYM:
          case MAT_STIFFNESS_SYM_L:
               return true;
          default:
               return false;
          }
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
     ElemBeam2(ElementTypes::TypeId eltype, octave_idx_type id, const Matrix& X, const Material* material, const int32NDArray& nodes, const BeamCrossSection& oSect, const ColumnVector& e2)
          :Element(eltype, id, X, material, nodes), R(3, 3) {

          FEM_ASSERT(X.rows() == 6);
          FEM_ASSERT(X.columns() == 2);
          FEM_ASSERT(nodes.numel() == 2);
          FEM_ASSERT(e2.numel() == 3);

          l = 0;

          for (octave_idx_type i = 0; i < 3; ++i) {
               double dX = X.xelem(i, 1) - X.xelem(i, 0);
               l += dX * dX;
               R.xelem(i, 0) = dX;
          }

          l = sqrt(l);

          if (l == 0) {
               throw std::runtime_error("zero beam length detected");
          }

          R.xelem(0, 2) = R.xelem(1, 0) * e2.xelem(2) - R.xelem(2, 0) * e2.xelem(1);
          R.xelem(1, 2) = R.xelem(2, 0) * e2.xelem(0) - R.xelem(0, 0) * e2.xelem(2);
          R.xelem(2, 2) = R.xelem(0, 0) * e2.xelem(1) - R.xelem(1, 0) * e2.xelem(0);

          R.xelem(0, 1) = R.xelem(1, 2) * R.xelem(2, 0) - R.xelem(2, 2) * R.xelem(1, 0);
          R.xelem(1, 1) = R.xelem(2, 2) * R.xelem(0, 0) - R.xelem(0, 2) * R.xelem(2, 0);
          R.xelem(2, 1) = R.xelem(0, 2) * R.xelem(1, 0) - R.xelem(1, 2) * R.xelem(0, 0);

          for (octave_idx_type j = 0; j < 3; ++j) {
               double n = 0;

               for (octave_idx_type i = 0; i < 3; ++i) {
                    double Rij = R.xelem(i, j);
                    n += Rij * Rij;
               }

               n = sqrt(n);

               if (n == 0.) {
                    throw std::runtime_error("orientation of beam cross section is not valid");
               }

               for (octave_idx_type i = 0; i < 3; ++i) {
                    R.xelem(i, j) /= n;
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

     virtual void Assemble(MatrixAss& mat, MeshInfo& info, const DofMap& dof, FemMatrixType eMatType) const {
          void (ElemBeam2::*pfn)(Matrix&) const;

          switch (eMatType) {
          case MAT_STIFFNESS:
          case MAT_STIFFNESS_SYM:
          case MAT_STIFFNESS_SYM_L:
               pfn = &ElemBeam2::StiffnessMatrix;
               break;

          case MAT_MASS:
          case MAT_MASS_SYM:
          case MAT_MASS_SYM_L:
               pfn = &ElemBeam2::MassMatrix;
               break;

          default:
               return;
          }

          Matrix A(12, 12, 0.), RA(12, 12);

          (this->*pfn)(A);

          for (octave_idx_type i0 = 0; i0 < 4; ++i0) {
               for (octave_idx_type j0 = 0; j0 < 4; ++j0) {
                    for (octave_idx_type i1 = 0; i1 < 3; ++i1) {
                         for (octave_idx_type j1 = 0; j1 < 3; ++j1) {
                              double RAij = 0;

                              for (octave_idx_type k = 0; k < 3; ++k) {
                                   RAij += R.xelem(i1, k) * A.xelem(3 * i0 + k, 3 * j0 + j1);
                              }

                              RA.xelem(3 * i0 + i1, 3 * j0 + j1) = RAij;
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
                                   Aij += R.xelem(j1, k) * RA.xelem(3 * i0 + i1, 3 * j0 + k);
                              }

                              A.xelem(3 * i0 + i1, 3 * j0 + j1) = Aij;
                         }
                    }
               }
          }

          for (octave_idx_type j = 0; j < 12; ++j) {
               for (octave_idx_type i = 0; i < j; ++i) {
                    A.xelem(j, i) = A.xelem(i, j);
               }
          }

          int32NDArray ndofidx(dim_vector(nodes.numel() * 6, 1), -1);

          for (octave_idx_type inode = 0; inode < nodes.numel(); ++inode) {
               for (octave_idx_type idof = 0; idof < 6; ++idof) {
                    ndofidx.xelem(inode * 6 + idof) = dof.GetNodeDofIndex(nodes.xelem(inode).value() - 1, idof);
               }
          }

          mat.Insert(A, ndofidx, ndofidx);
     }

     virtual octave_idx_type iGetWorkSpaceSize(FemMatrixType eMatType) const {
          return 12 * 12;
     }

     virtual double dGetMass() const {
          return rhoA * l;
     }

     void StiffnessMatrix(Matrix& Ke) const {
          const double l2 = l * l;
          const double l3 = l2 * l;

          Ke.xelem(0,0) = (EA)/l;
          Ke.xelem(1,1) = (12*EIz)/((12*ky+1)*l3);
          Ke.xelem(2,2) = (12*EIy)/((12*kz+1)*l3);
          Ke.xelem(3,3) = (GIt)/l;
          Ke.xelem(2,4) = -(6*EIy)/((12*kz+1)*l2);
          Ke.xelem(4,4) = (4*EIy*(3*kz+1))/((12*kz+1)*l);
          Ke.xelem(1,5) = (6*EIz)/((12*ky+1)*l2);
          Ke.xelem(5,5) = (4*EIz*(3*ky+1))/((12*ky+1)*l);
          Ke.xelem(0,6) = -(EA)/l;
          Ke.xelem(6,6) = (EA)/l;
          Ke.xelem(1,7) = -(12*EIz)/((12*ky+1)*l3);
          Ke.xelem(5,7) = -(6*EIz)/((12*ky+1)*l2);
          Ke.xelem(7,7) = (12*EIz)/((12*ky+1)*l3);
          Ke.xelem(2,8) = -(12*EIy)/((12*kz+1)*l3);
          Ke.xelem(4,8) = (6*EIy)/((12*kz+1)*l2);
          Ke.xelem(8,8) = (12*EIy)/((12*kz+1)*l3);
          Ke.xelem(3,9) = -(GIt)/l;
          Ke.xelem(9,9) = (GIt)/l;
          Ke.xelem(2,10) = -(6*EIy)/((12*kz+1)*l2);
          Ke.xelem(4,10) = -(2*EIy*(6*kz-1))/((12*kz+1)*l);
          Ke.xelem(8,10) = (6*EIy)/((12*kz+1)*l2);
          Ke.xelem(10,10) = (4*EIy*(3*kz+1))/((12*kz+1)*l);
          Ke.xelem(1,11) = (6*EIz)/((12*ky+1)*l2);
          Ke.xelem(5,11) = -(2*EIz*(6*ky-1))/((12*ky+1)*l);
          Ke.xelem(7,11) = -(6*EIz)/((12*ky+1)*l2);
          Ke.xelem(11,11) = (4*EIz*(3*ky+1))/((12*ky+1)*l);
     }

     void MassMatrix(Matrix& Me) const {
          const double l2 = l * l;
          const double ky2 = ky * ky;
          const double kz2 = kz * kz;
          const double b1 = (12*ky+1);
          const double b2 = (12*kz+1);
          const double a1 = b1 * b1;
          const double a2 = b2 * b2;

          Me.xelem(0,0) = (l*rhoA)/3.0E+0;
          Me.xelem(1,1) = ((42*rhoIz+1680*ky2*l2*rhoA+294*ky*l2*rhoA+13*l2*rhoA)/(a1*l))/3.5E+1;
          Me.xelem(2,2) = ((42*rhoIy+1680*kz2*l2*rhoA+294*kz*l2*rhoA+13*l2*rhoA)/(a2*l))/3.5E+1;
          Me.xelem(3,3) = (l*rhoIp)/3.0E+0;
          Me.xelem(2,4) = ((1260*kz*rhoIy-21*rhoIy-1260*kz2*l2*rhoA-231*kz*l2*rhoA-11*l2*rhoA)/a2)/2.1E+2;
          Me.xelem(4,4) = ((l*(5040*kz2*rhoIy+210*kz*rhoIy+14*rhoIy+126*kz2*l2*rhoA+21*kz*l2*rhoA+l2*rhoA))/a2)/1.05E+2;
          Me.xelem(1,5) = -((1260*ky*rhoIz-21*rhoIz-1260*ky2*l2*rhoA-231*ky*l2*rhoA-11*l2*rhoA)/a1)/2.1E+2;
          Me.xelem(5,5) = ((l*(5040*ky2*rhoIz+210*ky*rhoIz+14*rhoIz+126*ky2*l2*rhoA+21*ky*l2*rhoA+l2*rhoA))/a1)/1.05E+2;
          Me.xelem(0,6) = (l*rhoA)/6.0E+0;
          Me.xelem(6,6) = (l*rhoA)/3.0E+0;
          Me.xelem(1,7) = ((-3.0E+0)*(28*rhoIz-560*ky2*l2*rhoA-84*ky*l2*rhoA-3*l2*rhoA))/(7.0E+1*a1*l);
          Me.xelem(5,7) = ((2520*ky*rhoIz-42*rhoIz+2520*ky2*l2*rhoA+378*ky*l2*rhoA+13*l2*rhoA)/a1)/4.2E+2;
          Me.xelem(7,7) = ((42*rhoIz+1680*ky2*l2*rhoA+294*ky*l2*rhoA+13*l2*rhoA)/(a1*l))/3.5E+1;
          Me.xelem(2,8) = ((-3.0E+0)*(28*rhoIy-560*kz2*l2*rhoA-84*kz*l2*rhoA-3*l2*rhoA))/(7.0E+1*a2*l);
          Me.xelem(4,8) = -((2520*kz*rhoIy-42*rhoIy+2520*kz2*l2*rhoA+378*kz*l2*rhoA+13*l2*rhoA)/a2)/4.2E+2;
          Me.xelem(8,8) = ((42*rhoIy+1680*kz2*l2*rhoA+294*kz*l2*rhoA+13*l2*rhoA)/(a2*l))/3.5E+1;
          Me.xelem(3,9) = (l*rhoIp)/6.0E+0;
          Me.xelem(9,9) = (l*rhoIp)/3.0E+0;
          Me.xelem(2,10) = ((2520*kz*rhoIy-42*rhoIy+2520*kz2*l2*rhoA+378*kz*l2*rhoA+13*l2*rhoA)/a2)/4.2E+2;
          Me.xelem(4,10) = ((l*(10080*kz2*rhoIy-840*kz*rhoIy-14*rhoIy-504*kz2*l2*rhoA-84*kz*l2*rhoA-3*l2*rhoA))/a2)/4.2E+2;
          Me.xelem(8,10) = -((1260*kz*rhoIy-21*rhoIy-1260*kz2*l2*rhoA-231*kz*l2*rhoA-11*l2*rhoA)/a2)/2.1E+2;
          Me.xelem(10,10) = ((l*(5040*kz2*rhoIy+210*kz*rhoIy+14*rhoIy+126*kz2*l2*rhoA+21*kz*l2*rhoA+l2*rhoA))/a2)/1.05E+2;
          Me.xelem(1,11) = -((2520*ky*rhoIz-42*rhoIz+2520*ky2*l2*rhoA+378*ky*l2*rhoA+13*l2*rhoA)/a1)/4.2E+2;
          Me.xelem(5,11) = ((l*(10080*ky2*rhoIz-840*ky*rhoIz-14*rhoIz-504*ky2*l2*rhoA-84*ky*l2*rhoA-3*l2*rhoA))/a1)/4.2E+2;
          Me.xelem(7,11) = ((1260*ky*rhoIz-21*rhoIz-1260*ky2*l2*rhoA-231*ky*l2*rhoA-11*l2*rhoA)/a1)/2.1E+2;
          Me.xelem(11,11) = ((l*(5040*ky2*rhoIz+210*ky*rhoIz+14*rhoIz+126*ky2*l2*rhoA+21*ky*l2*rhoA+l2*rhoA))/a1)/1.05E+2;
     }

private:
     Matrix R;
     double l, ky, kz, EA, GAy, GAz, GIt, EIy, EIz, rhoA, rhoIp, rhoIy, rhoIz;
};

class Element3D: public Element
{
public:
     Element3D(ElementTypes::TypeId eltype, octave_idx_type id, const Matrix& X, const Material* material, const int32NDArray& nodes, const StrainField& oRefStrain)
          :Element(eltype, id, X, material, nodes) {

          FEM_ASSERT(X.rows() == 3);

          const octave_idx_type iNumNodes = nodes.numel();
          const octave_idx_type iNumLoadsTemp = oRefStrain.rgTemperature.numel();
          const octave_idx_type iNumLoadsStrain = oRefStrain.rgRefStrain.numel();
          
          iNumPreLoads = std::max(iNumLoadsTemp, iNumLoadsStrain);
          
          FEM_ASSERT(iNumLoadsTemp && iNumLoadsStrain ? iNumLoadsTemp == iNumLoadsStrain : true);
          
          if (iNumLoadsTemp) {
               dTheta.resize(iNumNodes, iNumLoadsTemp);
               
               for (octave_idx_type j = 0; j < iNumLoadsTemp; ++j) {
                    const NDArray dThetaj = oRefStrain.rgTemperature.xelem(j).array_value();

                    for (octave_idx_type i = 0; i < iNumNodes; ++i) {
                         dTheta.xelem(i, j) = dThetaj.xelem(nodes.xelem(i).value() - 1);
                    }
               }
          }

          if (iNumLoadsStrain) {
               const octave_idx_type iNumStrains = material->LinearElasticity().rows();
               
               epsilonRef.resize(dim_vector(iNumStrains, iNumNodes, iNumLoadsStrain), 0.);

               for (octave_idx_type k = 0; k < iNumLoadsStrain; ++k) {
                    const octave_scalar_map maEpsilonRefk = oRefStrain.rgRefStrain.xelem(k).scalar_map_value();

                    const std::string strElemName = ElementTypes::GetType(eltype).name;
                    const auto iterEpsilonRefk = maEpsilonRefk.seek(strElemName);

                    if (iterEpsilonRefk == maEpsilonRefk.end()) {
                         continue;
                    }

                    const NDArray epsilonRefk = maEpsilonRefk.contents(iterEpsilonRefk).array_value();
                    
                    if (epsilonRefk.ndims() != 3 || epsilonRefk.dim1() < id || epsilonRefk.dim2() != iNumNodes || epsilonRefk.dim3() != iNumStrains) {
                         throw std::runtime_error("fem_ass_matrix: invalid number of dimensions for load_case.epsilon0."s + strElemName);
                    }
                    
                    for (octave_idx_type j = 0; j < iNumNodes; ++j) {
                         for (octave_idx_type i = 0; i < iNumStrains; ++i) {
                              epsilonRef.xelem(i, j, k) = epsilonRefk.xelem(id - 1, j, i);
                         }
                    }
               }
          }
     }

     virtual void Assemble(MatrixAss& mat, MeshInfo& info, const DofMap& dof, FemMatrixType eMatType) const {
          void (Element3D::*pFunc)(Matrix&, MeshInfo&, FemMatrixType) const;

          const octave_idx_type iNumDof = iGetNumDof(eMatType);

          octave_idx_type iNumRows, iNumCols;

          switch (eMatType) {
          case MAT_STIFFNESS:
          case MAT_STIFFNESS_SYM:
          case MAT_STIFFNESS_SYM_L:
               pFunc = &Element3D::StiffnessMatrix;
               iNumRows = iNumCols = iNumDof;
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
               pFunc = &Element3D::ThermalLoadVector;
               iNumRows = iNumDof;
               iNumCols = iNumPreLoads;
               break;

          case MAT_ACCEL_LOAD:
               pFunc = &Element3D::AccelerationLoad;
               iNumRows = iNumDof;
               iNumCols = 3;
               break;
               
          case MAT_THERMAL_COND:
               pFunc = &Element3D::ThermalConductivityMatrix;
               iNumRows = iNumCols = iNumDof;
               break;

          case MAT_HEAT_CAPACITY:
               pFunc = &Element3D::HeatCapacityMatrix;
               iNumRows = iNumCols = iNumDof;
               break;
               
          default:
               return;
          }

          if (iNumCols == 0) {
               return;
          }

          int32NDArray dofidx(dim_vector(iNumDof, 1), 0);

          const octave_idx_type inodemaxdof = (eMatType & DofMap::DO_THERMAL) ? 1 : 3;
          
          for (octave_idx_type inode = 0; inode < nodes.numel(); ++inode) {
               for (octave_idx_type idof = 0; idof < inodemaxdof; ++idof) {
                    dofidx.xelem(inode * inodemaxdof + idof) = dof.GetNodeDofIndex(nodes.xelem(inode).value() - 1, idof);
               }
          }

          Matrix Ae(iNumRows, iNumCols, 0.);

          (this->*pFunc)(Ae, info, eMatType);

          switch (eMatType) {
          case MAT_ACCEL_LOAD:
          case VEC_LOAD_CONSISTENT:
          case VEC_LOAD_LUMPED: {
               int32NDArray dofidxcol(dim_vector(iNumCols, 1));

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
          switch (eMatType) {
          case MAT_THERMAL_COND:
          case MAT_HEAT_CAPACITY:
               return nodes.numel();
          default:
               return nodes.numel() * 3;
          }
     }

     virtual octave_idx_type iGetWorkSpaceSize(FemMatrixType eMatType) const {
          switch (eMatType) {
          case MAT_MASS:
          case MAT_STIFFNESS:
          case MAT_DAMPING:
          case MAT_THERMAL_COND:
          case MAT_HEAT_CAPACITY: {
               const octave_idx_type iNumDof = iGetNumDof(eMatType);

               return iNumDof * iNumDof;
          }
          case MAT_MASS_SYM:
          case MAT_MASS_SYM_L:
          case MAT_STIFFNESS_SYM:
          case MAT_STIFFNESS_SYM_L:
          case MAT_DAMPING_SYM:
          case MAT_DAMPING_SYM_L: {
               const octave_idx_type iNumDof = iGetNumDof(eMatType);

               return (iNumDof + 1) * iNumDof / 2;
          }
          case MAT_MASS_LUMPED:
               return iGetNumDof(eMatType);

          case MAT_ACCEL_LOAD:
               return iGetNumDof(eMatType) * 3;

          default:
               return 0;
          }
     }

     virtual void PostProcElem(NDArray& mat, FemMatrixType eMatType, const NDArray& U) const {
          switch (eMatType) {
          case VEC_INERTIA_M1:
               InertiaMoment1(mat, eMatType);
               break;

          case MAT_INERTIA_J:
               InertiaMatrix(mat, eMatType);
               break;

          case MAT_INERTIA_INV3:
               InertiaInv3(mat, eMatType, U);
               break;

          case MAT_INERTIA_INV4:
               InertiaInv4(mat, eMatType, U);
               break;

          case MAT_INERTIA_INV5:
               InertiaInv5(mat, eMatType, U);
               break;

          case MAT_INERTIA_INV8:
               InertiaInv8(mat, eMatType, U);
               break;

          case MAT_INERTIA_INV9:
               InertiaInv9(mat, eMatType, U);
               break;

          case VEC_STRESS_CAUCH:
               StressNodalElem(mat, eMatType, U);
               break;

          case VEC_STRAIN_TOTAL:
               StrainNodalElem(mat, eMatType, U);
               break;
               
          default:
               break;
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

     virtual double dGetMass() const {
          return dGetVolume() * material->Density();
     }

protected:
     virtual const IntegrationRule& GetIntegrationRule(FemMatrixType eMatType) const=0;
     virtual double Jacobian(const ColumnVector& rv, Matrix& J) const=0;
     virtual void StrainMatrix(const ColumnVector& rv, const Matrix& invJ, Matrix& B) const=0;
     virtual void DispInterpMatrix(const ColumnVector& rv, Matrix& H) const=0;
     virtual Matrix InterpGaussToNodal(FemMatrixType eMatType, const Matrix& taun) const=0;
     virtual void ScalarInterpMatrix(const ColumnVector& rv, Matrix& Hs, octave_idx_type irow) const=0;
     virtual void TemperatureGradientMatrix(const ColumnVector& rv, const Matrix& invJ, Matrix& Bt) const=0;
     
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

          Matrix J(iNumDir, iNumDir), B(iNumStrains, iNumDof), CB(iNumStrains, iNumDof);

          for (octave_idx_type i = 0; i < oIntegRule.iGetNumEvalPoints(); ++i) {
               const double alpha = oIntegRule.dGetWeight(i);

               for (octave_idx_type j = 0; j < iNumDir; ++j) {
                    rv.xelem(j) = oIntegRule.dGetPosition(i, j);
               }

               const double detJ = Jacobian(rv, J);

               AddMeshInfo(info, oIntegRule, detJ);

               StrainMatrix(rv, J.inverse(), B);

               for (octave_idx_type l = 0; l < iNumStrains; ++l) {
                    for (octave_idx_type m = 0; m < iNumDof; ++m) {
                         double CBlm = 0.;

                         for (octave_idx_type n = 0; n < iNumStrains; ++n) {
                              CBlm += detJ * alpha * C.xelem(l, n) * B.xelem(n, m);
                         }

                         FEM_ASSERT(std::isfinite(CBlm));

                         CB.xelem(l, m) = CBlm;
                    }
               }

#ifdef DEBUG
               for (octave_idx_type i = 0; i < CB.rows(); ++i) {
                    for (octave_idx_type j = 0; j < CB.columns(); ++j) {
                         FEM_ASSERT(std::isfinite(CB(i, j)));
                    }
               }
#endif

               for (octave_idx_type l = 0; l < iNumDof; ++l) {
                    for (octave_idx_type m = l; m < iNumDof; ++m) {
                         double Kelm = 0.;

                         for (octave_idx_type n = 0; n < iNumStrains; ++n) {
                              Kelm += B.xelem(n, l) * CB.xelem(n, m);
                         }

                         FEM_ASSERT(std::isfinite(Kelm));

                         Ke.xelem(l, m) += Kelm;
                    }
               }
          }

          for (octave_idx_type i = 1; i < iNumDof; ++i) {
               for (octave_idx_type j = 0; j < i; ++j) {
                    Ke.xelem(i, j) = Ke.xelem(j, i);
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
                              Melm += H.xelem(n, l) * H.xelem(n, m);
                         }

                         FEM_ASSERT(std::isfinite(Melm));

                         Me.xelem(l, m) += Melm * alpha * rho * detJ;
                    }
               }
          }

          for (octave_idx_type i = 1; i < iNumDof; ++i) {
               for (octave_idx_type j = 0; j < i; ++j) {
                    Me.xelem(i, j) = Me.xelem(j, i);
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
          const double alpha = material->AlphaDamping();

          if (alpha) {
               MassMatrix(De, info, static_cast<FemMatrixType>(MAT_MASS | (eMatType & MAT_SYM_MASK)));

               for (octave_idx_type j = 0; j < De.columns(); ++j) {
                    for (octave_idx_type i = 0; i < De.rows(); ++i) {
                         De.xelem(i, j) *= alpha;
                    }
               }
          }

          const double beta = material->BetaDamping();

          if (beta) {
               Matrix Ke(De.rows(), De.columns(), 0.);

               StiffnessMatrix(Ke, info, static_cast<FemMatrixType>(MAT_STIFFNESS | (eMatType & MAT_SYM_MASK)));

               for (octave_idx_type j = 0; j < De.columns(); ++j) {
                    for (octave_idx_type i = 0; i < De.rows(); ++i) {
                         De.xelem(i, j) += beta * Ke.xelem(i, j);
                    }
               }
          }
     }

     void AccelerationLoad(Matrix& C1, MeshInfo& info, FemMatrixType eMatType) const {
          const IntegrationRule& oIntegRule = GetIntegrationRule(eMatType);
          const octave_idx_type iNumDof = iGetNumDof(eMatType);
          const octave_idx_type iNumDir = oIntegRule.iGetNumDirections();
          ColumnVector rv(iNumDir);
          const double rho = material->Density();
          const octave_idx_type iNumDisp = X.rows();

          FEM_ASSERT(C1.rows() == iNumDof);
          FEM_ASSERT(C1.columns() == iNumDisp);

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

               for (octave_idx_type k = 0; k < iNumDisp; ++k) {
                    for (octave_idx_type j = 0; j < iNumDof; ++j) {
                         C1.xelem(j, k) += H.xelem(k, j) * dm;
                    }
               }
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

          FEM_ASSERT(S.ndims() == 2);
          FEM_ASSERT(S.rows() == 3);
          FEM_ASSERT(S.columns() == 1);
          FEM_ASSERT(iNumDisp == S.rows());

          Matrix J(iNumDir, iNumDir), H(iNumDisp, iNumDof);

          for (octave_idx_type i = 0; i < oIntegRule.iGetNumEvalPoints(); ++i) {
               const double alpha = oIntegRule.dGetWeight(i);

               for (octave_idx_type j = 0; j < iNumDir; ++j) {
                    rv.xelem(j) = oIntegRule.dGetPosition(i, j);
               }

               const double detJ = Jacobian(rv, J);;

               DispInterpMatrix(rv, H);

               for (octave_idx_type l = 0; l < S.rows(); ++l) {
                    double fil = 0.;

                    for (octave_idx_type n = 0; n < iNumNodes; ++n) {
                         for (octave_idx_type m = 0; m < iNumDisp; ++m) {
                              fil += H.xelem(l, iNumDisp * n + m) * X.xelem(m, n);
                         }
                    }

                    S.xelem(l) += fil * alpha * rho * detJ;
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
                              fi.xelem(l) += H.xelem(l, iNumDisp * n + m) * X.xelem(m, n);
                         }
                    }
               }

               const double dmi = alpha * rho * detJ;

               Inv7.xelem(0, 0) += (fi.xelem(1) * fi.xelem(1) + fi.xelem(2) * fi.xelem(2)) * dmi;
               Inv7.xelem(0, 1) -= (fi.xelem(0) * fi.xelem(1)) * dmi;
               Inv7.xelem(0, 2) -= (fi.xelem(0) * fi.xelem(2)) * dmi;
               Inv7.xelem(1, 1) += (fi.xelem(0) * fi.xelem(0) + fi.xelem(2) * fi.xelem(2)) * dmi;
               Inv7.xelem(1, 2) -= (fi.xelem(1) * fi.xelem(2)) * dmi;
               Inv7.xelem(2, 2) += (fi.xelem(0) * fi.xelem(0) + fi.xelem(1) * fi.xelem(1)) * dmi;
          }
     }

     void InertiaInv3(NDArray& Inv3, FemMatrixType eMatType, const NDArray& U) const {
          const IntegrationRule& oIntegRule = GetIntegrationRule(eMatType);
          const octave_idx_type iNumDof = iGetNumDof(eMatType);
          const octave_idx_type iNumDir = oIntegRule.iGetNumDirections();
          const double rho = material->Density();
          const octave_idx_type iNumDisp = X.rows();
          const octave_idx_type iNumNodes = nodes.numel();
          ColumnVector rv(iNumDir);

          FEM_ASSERT(U.ndims() == 3);
          FEM_ASSERT(U.dim2() >= iNumDisp);
          FEM_ASSERT(Inv3.ndims() == 2);
          FEM_ASSERT(Inv3.rows() == 3);
          FEM_ASSERT(Inv3.columns() == U.dim3());
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

               for (octave_idx_type j = 0; j < U.dim3(); ++j) {
                    for (octave_idx_type l = 0; l < iNumDisp; ++l) {
                         double Uil = 0.;

                         for (octave_idx_type m = 0; m < iNumDisp; ++m) {
                              for (octave_idx_type n = 0; n < iNumNodes; ++n) {
                                   Uil += H.xelem(l, iNumDisp * n + m) * U.xelem(nodes.xelem(n).value() - 1, m, j);
                              }
                         }

                         Inv3.xelem(l, j) += dmi * Uil;
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
          ColumnVector rv(iNumDir);

          FEM_ASSERT(U.ndims() == 3);
          FEM_ASSERT(U.dim2() >= iNumDisp);
          FEM_ASSERT(Inv4.rows() == 3);
          FEM_ASSERT(Inv4.columns() == U.dim3());
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

               for (octave_idx_type j = 0; j < U.dim3(); ++j) {
                    for (octave_idx_type l = 0; l < iNumDisp; ++l) {
                         Ui.xelem(l) = 0.;
                         fi.xelem(l) = 0.;

                         for (octave_idx_type m = 0; m < iNumDisp; ++m) {
                              for (octave_idx_type n = 0; n < iNumNodes; ++n) {
                                   const double Hlmn = H.xelem(l, iNumDisp * n + m);
                                   Ui.xelem(l) += Hlmn * U.xelem(nodes.xelem(n).value() - 1, m, j);
                                   fi.xelem(l) += Hlmn * X.xelem(m, n);
                              }
                         }
                    }

                    Inv4.xelem(0, j) += (fi.xelem(1) * Ui.xelem(2) - Ui.xelem(1) * fi.xelem(2)) * dmi;
                    Inv4.xelem(1, j) += (Ui.xelem(0) * fi.xelem(2) - fi.xelem(0) * Ui.xelem(2)) * dmi;
                    Inv4.xelem(2, j) += (fi.xelem(0) * Ui.xelem(1) - Ui.xelem(0) * fi.xelem(1)) * dmi;
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
          ColumnVector rv(iNumDir);

          FEM_ASSERT(U.ndims() == 3);
          FEM_ASSERT(U.dim2() >= iNumDisp);
          FEM_ASSERT(Inv5.ndims() == 3);
          FEM_ASSERT(Inv5.dim1() == 3);
          FEM_ASSERT(Inv5.dim2() == U.dim3());
          FEM_ASSERT(Inv5.dim3() == U.dim3());
          FEM_ASSERT(iNumDisp == Inv5.dim1());

          Matrix J(iNumDir, iNumDir), H(iNumDisp, iNumDof);
          NDArray Ui(dim_vector(iNumDisp, U.dim3(), iNumGauss), 0.);
          ColumnVector dmi(iNumGauss);

          for (octave_idx_type i = 0; i < iNumGauss; ++i) {
               const double alpha = oIntegRule.dGetWeight(i);

               for (octave_idx_type j = 0; j < iNumDir; ++j) {
                    rv.xelem(j) = oIntegRule.dGetPosition(i, j);
               }

               const double detJ = Jacobian(rv, J);

               DispInterpMatrix(rv, H);

               dmi.xelem(i) = alpha * rho * detJ;

               for (octave_idx_type j = 0; j < U.dim3(); ++j) {
                    for (octave_idx_type l = 0; l < iNumDisp; ++l) {
                         for (octave_idx_type m = 0; m < iNumDisp; ++m) {
                              for (octave_idx_type n = 0; n < iNumNodes; ++n) {
                                   Ui.xelem(l, j, i) += H.xelem(l, iNumDisp * n + m) * U.xelem(nodes.xelem(n).value() - 1, m, j);
                              }
                         }
                    }
               }
          }

          for (octave_idx_type i = 0; i < iNumGauss; ++i) {
               for (octave_idx_type j = 0; j < U.dim3(); ++j) {
                    for (octave_idx_type k = 0; k < U.dim3(); ++k) {
                         Inv5.xelem(0, k, j) += dmi.xelem(i) * (Ui.xelem(1, j, i) * Ui.xelem(2, k, i) - Ui.xelem(2, j, i) * Ui.xelem(1, k, i));
                         Inv5.xelem(1, k, j) += dmi.xelem(i) * (Ui.xelem(2, j, i) * Ui.xelem(0, k, i) - Ui.xelem(0, j, i) * Ui.xelem(2, k, i));
                         Inv5.xelem(2, k, j) += dmi.xelem(i) * (Ui.xelem(0, j, i) * Ui.xelem(1, k, i) - Ui.xelem(1, j, i) * Ui.xelem(0, k, i));
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
          ColumnVector rv(iNumDir);

          FEM_ASSERT(U.ndims() == 3);
          FEM_ASSERT(U.dim2() >= iNumDisp);
          FEM_ASSERT(Inv8.ndims() == 3);
          FEM_ASSERT(Inv8.dim1() == 3);
          FEM_ASSERT(Inv8.dim2() == 3);
          FEM_ASSERT(Inv8.dim3() == U.dim3());
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

               for (octave_idx_type j = 0; j < U.dim3(); ++j) {
                    for (octave_idx_type l = 0; l < iNumDisp; ++l) {
                         Ui.xelem(l) = 0.;
                         fi.xelem(l) = 0.;

                         for (octave_idx_type m = 0; m < iNumDisp; ++m) {
                              for (octave_idx_type n = 0; n < iNumNodes; ++n) {
                                   const double Hlmn = H.xelem(l, iNumDisp * n + m);
                                   Ui.xelem(l) += Hlmn * U.xelem(nodes.xelem(n).value() - 1, m, j);
                                   fi.xelem(l) += Hlmn * X.xelem(m, n);
                              }
                         }
                    }

                    Inv8.xelem(0,0,j) += (fi.xelem(2)*Ui.xelem(2)+fi.xelem(1)*Ui.xelem(1))*dmi;
                    Inv8.xelem(1,0,j) += -fi.xelem(0)*Ui.xelem(1)*dmi;
                    Inv8.xelem(2,0,j) += -fi.xelem(0)*Ui.xelem(2)*dmi;
                    Inv8.xelem(0,1,j) += -Ui.xelem(0)*fi.xelem(1)*dmi;
                    Inv8.xelem(1,1,j) += (fi.xelem(2)*Ui.xelem(2)+fi.xelem(0)*Ui.xelem(0))*dmi;
                    Inv8.xelem(2,1,j) += -fi.xelem(1)*Ui.xelem(2)*dmi;
                    Inv8.xelem(0,2,j) += -Ui.xelem(0)*fi.xelem(2)*dmi;
                    Inv8.xelem(1,2,j) += -Ui.xelem(1)*fi.xelem(2)*dmi;
                    Inv8.xelem(2,2,j) += (fi.xelem(1)*Ui.xelem(1)+fi.xelem(0)*Ui.xelem(0))*dmi;
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
          const octave_idx_type iNumGauss = oIntegRule.iGetNumEvalPoints();
          ColumnVector rv(iNumDir);

          FEM_ASSERT(U.ndims() == 3);
          FEM_ASSERT(U.dim2() >= iNumDisp);
          FEM_ASSERT(Inv9.ndims() == 4);
          FEM_ASSERT(Inv9.dim1() == 3);
          FEM_ASSERT(Inv9.dim2() == 3);
          FEM_ASSERT(Inv9.dim3() == U.dim3());
          FEM_ASSERT(Inv9.dims()(3) == U.dim3());
          FEM_ASSERT(iNumDisp == Inv9.dim1());

          Matrix J(iNumDir, iNumDir), H(iNumDisp, iNumDof);
          NDArray Ui(dim_vector(iNumDisp, U.dim3(), iNumGauss), 0.);
          ColumnVector dmi(iNumGauss);

          for (octave_idx_type i = 0; i < iNumGauss; ++i) {
               const double alpha = oIntegRule.dGetWeight(i);

               for (octave_idx_type j = 0; j < iNumDir; ++j) {
                    rv.xelem(j) = oIntegRule.dGetPosition(i, j);
               }

               const double detJ = Jacobian(rv, J);

               DispInterpMatrix(rv, H);

               dmi.xelem(i) = alpha * rho * detJ;

               for (octave_idx_type j = 0; j < U.dim3(); ++j) {
                    for (octave_idx_type l = 0; l < iNumDisp; ++l) {
                         for (octave_idx_type m = 0; m < iNumDisp; ++m) {
                              for (octave_idx_type n = 0; n < iNumNodes; ++n) {
                                   Ui.xelem(l, j, i) += H.xelem(l, iNumDisp * n + m) * U.xelem(nodes.xelem(n).value() - 1, m, j);
                              }
                         }
                    }
               }
          }

          for (octave_idx_type i = 0; i < iNumGauss; ++i) {
               for (octave_idx_type j = 0; j < U.dim3(); ++j) {
                    for (octave_idx_type k = 0; k < U.dim3(); ++k) {
                         const octave_idx_type jk = j + k * Inv9.dim3();

                         Inv9.xelem(0, 0, jk) += (-Ui.xelem(2,j,i)*Ui.xelem(2,k,i)-Ui.xelem(1,j,i)*Ui.xelem(1,k,i))*dmi.xelem(i);
                         Inv9.xelem(1, 0, jk) += Ui.xelem(0,j,i)*Ui.xelem(1,k,i)*dmi.xelem(i);
                         Inv9.xelem(2, 0, jk) += Ui.xelem(0,j,i)*Ui.xelem(2,k,i)*dmi.xelem(i);
                         Inv9.xelem(0, 1, jk) += Ui.xelem(0,k,i)*Ui.xelem(1,j,i)*dmi.xelem(i);
                         Inv9.xelem(1, 1, jk) += (-Ui.xelem(2,j,i)*Ui.xelem(2,k,i)-Ui.xelem(0,j,i)*Ui.xelem(0,k,i))*dmi.xelem(i);
                         Inv9.xelem(2, 1, jk) += Ui.xelem(1,j,i)*Ui.xelem(2,k,i)*dmi.xelem(i);
                         Inv9.xelem(0, 2, jk) += Ui.xelem(0,k,i)*Ui.xelem(2,j,i)*dmi.xelem(i);
                         Inv9.xelem(1, 2, jk) += Ui.xelem(1,k,i)*Ui.xelem(2,j,i)*dmi.xelem(i);
                         Inv9.xelem(2, 2, jk) += (-Ui.xelem(1,j,i)*Ui.xelem(1,k,i)-Ui.xelem(0,j,i)*Ui.xelem(0,k,i))*dmi.xelem(i);
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
          FEM_ASSERT(iNumDof == 3 * X.columns());

          Matrix J(iNumDir, iNumDir), B(iNumStrains, iNumDof);
          ColumnVector Ue(iNumDof);
          Matrix epsilong(iNumGauss, iNumStrains * iNumLoads);

          for (octave_idx_type i = 0; i < iNumGauss; ++i) {
               for (octave_idx_type j = 0; j < iNumDir; ++j) {
                    rv.xelem(j) = oIntegRule.dGetPosition(i, j);
               }

               Jacobian(rv, J);

               StrainMatrix(rv, J.inverse(), B);

               for (octave_idx_type l = 0; l < iNumLoads; ++l) {
                    for (octave_idx_type j = 0; j < iNumNodes; ++j) {
                         for (octave_idx_type k = 0; k < 3; ++k) {
                              FEM_ASSERT(nodes(j).value() > 0);
                              FEM_ASSERT(nodes(j).value() <= U.dim1());
                              Ue.xelem(3 * j + k) = U.xelem(nodes.xelem(j).value() - 1, k, l);
                         }
                    }
                    
                    for (octave_idx_type j = 0; j < iNumStrains; ++j) {
                         double epsilongj = 0.;

                         for (octave_idx_type k = 0; k < iNumDof; ++k) {
                              epsilongj += B.xelem(j, k) * Ue.xelem(k);
                         }
                         
                         epsilong.xelem(i, l * iNumStrains + j) = epsilongj;
                    }
               }
          }

          const Matrix epsilonen = InterpGaussToNodal(eMatType, epsilong);

          FEM_ASSERT(epsilonen.rows() == iNumNodes);
          FEM_ASSERT(epsilonen.columns() == epsilong.columns());

          for (octave_idx_type k = 0; k < iNumLoads; ++k) {
               for (octave_idx_type j = 0; j < iNumStrains; ++j) {
                    for (octave_idx_type i = 0; i < iNumNodes; ++i) {
                         epsilonn.xelem(id - 1, i, j + k * iNumStrains) = epsilonen.xelem(i, k * iNumStrains + j);
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
          FEM_ASSERT(iNumDof == 3 * X.columns());

          Matrix J(iNumDir, iNumDir), B(iNumStrains, iNumDof), CB(iNumStrains, iNumDof), Ht;
          ColumnVector Ue(iNumDof), epsilonik(iNumPreLoads ? iNumStrains : 0);
          Matrix taug(iNumGauss, iNumStrains * iNumLoads);

          if (iNumPreLoads) {
               Ht.resize(1, iNumNodes);
          }

          for (octave_idx_type i = 0; i < iNumGauss; ++i) {
               for (octave_idx_type j = 0; j < iNumDir; ++j) {
                    rv.xelem(j) = oIntegRule.dGetPosition(i, j);
               }

               Jacobian(rv, J);

               StrainMatrix(rv, J.inverse(), B);

               if (iNumPreLoads) {
                    ScalarInterpMatrix(rv, Ht, 0);
               }

               for (octave_idx_type j = 0; j < iNumStrains; ++j) {
                    for (octave_idx_type k = 0; k < iNumDof; ++k) {
                         double CBjk = 0;

                         for (octave_idx_type l = 0; l < iNumStrains; ++l) {
                              CBjk += C.xelem(j, l) * B.xelem(l, k);
                         }

                         CB.xelem(j, k) = CBjk;
                    }
               }

               for (octave_idx_type l = 0; l < iNumLoads; ++l) {
                    for (octave_idx_type j = 0; j < iNumNodes; ++j) {
                         for (octave_idx_type k = 0; k < 3; ++k) {
                              FEM_ASSERT(nodes(j).value() > 0);
                              FEM_ASSERT(nodes(j).value() <= U.dim1());
                              Ue.xelem(3 * j + k) = U.xelem(nodes.xelem(j).value() - 1, k, l);
                         }
                    }

                    double dThetail = 0.;

                    if (dTheta.columns()) {
                         for (octave_idx_type j = 0; j < iNumNodes; ++j) {
                              dThetail += Ht.xelem(j) * dTheta.xelem(j, l);
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
                                   epsilonik.xelem(k) += Ht.xelem(j) * epsilonRef.xelem(k, j, l);
                              }
                         }
                    }
                    
                    for (octave_idx_type j = 0; j < iNumStrains; ++j) {
                         double taugj = 0.;

                         for (octave_idx_type k = 0; k < iNumDof; ++k) {
                              taugj += CB.xelem(j, k) * Ue.xelem(k);
                         }

                         if (iNumPreLoads) {
                              for (octave_idx_type k = 0; k < iNumStrains; ++k) {
                                   taugj -= C.xelem(j, k) * epsilonik.xelem(k);
                              }
                         }
                         
                         taug.xelem(i, l * iNumStrains + j) = taugj;
                    }
               }
          }

          const Matrix tauen = InterpGaussToNodal(eMatType, taug);

          FEM_ASSERT(tauen.rows() == iNumNodes);
          FEM_ASSERT(tauen.columns() == taug.columns());

          for (octave_idx_type k = 0; k < iNumLoads; ++k) {
               for (octave_idx_type j = 0; j < iNumStrains; ++j) {
                    for (octave_idx_type i = 0; i < iNumNodes; ++i) {
                         taun.xelem(id - 1, i, j + k * iNumStrains) = tauen.xelem(i, k * iNumStrains + j);
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

          Matrix Ht(1, iNumNodes);
          ColumnVector rv(iNumDir);

          Matrix J(iNumDir, iNumDir), B(iNumStrains, iNumDof);
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
                              dThetaik += Ht.xelem(j) * dTheta.xelem(j, k);
                         }
                    }

                    FEM_ASSERT(iNumStrains >= 3);
                    
                    const double epsilonk = gamma * dThetaik;
                    
                    for (octave_idx_type i = 0; i < iNumStrains; ++i) {
                         epsilon.xelem(i, k) = i < 3 ? epsilonk : 0.;
                    }
                    
                    if (epsilonRef.numel()) {                         
                         for (octave_idx_type j = 0; j < iNumNodes; ++j) {
                              for (octave_idx_type i = 0; i < iNumStrains; ++i) {
                                   epsilon.xelem(i, k) += Ht.xelem(j) * epsilonRef.xelem(i, j, k);
                              }                             
                         }
                    }
               }

               for (octave_idx_type l = 0; l < iNumPreLoads; ++l) {
                    for (octave_idx_type j = 0; j < iNumStrains; ++j) {
                         double Cepsilonjl = 0.;

                         for (octave_idx_type k = 0; k < iNumStrains; ++k) {
                              Cepsilonjl += C.xelem(j, k) * epsilon.xelem(k, l);
                         }

                         Cepsilon.xelem(j, l) = Cepsilonjl;
                    }
               }

               const double detJ = Jacobian(rv, J);

               AddMeshInfo(info, oIntegRule, detJ);

               StrainMatrix(rv, J.inverse(), B);

               for (octave_idx_type j = 0; j < iNumPreLoads; ++j) {
                    for (octave_idx_type k = 0; k < iNumDof; ++k) {
                         double Rkj = 0.;

                         for (octave_idx_type l = 0; l < iNumStrains; ++l) {
                              Rkj += B.xelem(l, k) * Cepsilon.xelem(l, j);
                         }

                         R.xelem(k, j) += detJ * alpha * Rkj;
                    }
               }
          }
     }

     void ThermalConductivityMatrix(Matrix& Ke, MeshInfo& info, FemMatrixType eMatType) const {
          const IntegrationRule& oIntegRule = GetIntegrationRule(eMatType);
          const octave_idx_type iNumDof = iGetNumDof(eMatType);
          const octave_idx_type iNumDir = oIntegRule.iGetNumDirections();
          ColumnVector rv(iNumDir);
          const Matrix& k = material->ThermalConductivity();
          const octave_idx_type iNumStrains = k.rows();

          FEM_ASSERT(k.rows() == k.columns());

          Matrix J(iNumDir, iNumDir), B(iNumStrains, iNumDof), kB(iNumStrains, iNumDof);

          for (octave_idx_type i = 0; i < oIntegRule.iGetNumEvalPoints(); ++i) {
               const double alpha = oIntegRule.dGetWeight(i);

               for (octave_idx_type j = 0; j < iNumDir; ++j) {
                    rv.xelem(j) = oIntegRule.dGetPosition(i, j);
               }

               const double detJ = Jacobian(rv, J);

               AddMeshInfo(info, oIntegRule, detJ);

               TemperatureGradientMatrix(rv, J.inverse(), B);

               for (octave_idx_type l = 0; l < iNumStrains; ++l) {
                    for (octave_idx_type m = 0; m < iNumDof; ++m) {
                         double kBlm = 0.;

                         for (octave_idx_type n = 0; n < iNumStrains; ++n) {
                              kBlm += detJ * alpha * k.xelem(l, n) * B.xelem(n, m);
                         }

                         FEM_ASSERT(std::isfinite(kBlm));

                         kB.xelem(l, m) = kBlm;
                    }
               }

#ifdef DEBUG
               for (octave_idx_type i = 0; i < kB.rows(); ++i) {
                    for (octave_idx_type j = 0; j < kB.columns(); ++j) {
                         FEM_ASSERT(std::isfinite(kB(i, j)));
                    }
               }
#endif

               for (octave_idx_type l = 0; l < iNumDof; ++l) {
                    for (octave_idx_type m = l; m < iNumDof; ++m) {
                         double Kelm = 0.;

                         for (octave_idx_type n = 0; n < iNumStrains; ++n) {
                              Kelm += B.xelem(n, l) * kB.xelem(n, m);
                         }

                         FEM_ASSERT(std::isfinite(Kelm));

                         Ke.xelem(l, m) += Kelm;
                    }
               }
          }

          for (octave_idx_type i = 1; i < iNumDof; ++i) {
               for (octave_idx_type j = 0; j < i; ++j) {
                    Ke.xelem(i, j) = Ke.xelem(j, i);
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

     void HeatCapacityMatrix(Matrix& Ce, MeshInfo& info, FemMatrixType eMatType) const {
          const IntegrationRule& oIntegRule = GetIntegrationRule(eMatType);
          const octave_idx_type iNumDof = iGetNumDof(eMatType);
          const octave_idx_type iNumDir = oIntegRule.iGetNumDirections();
          ColumnVector rv(iNumDir);
          const double rhocp = material->Density() * material->HeatCapacity();

          Matrix J(iNumDir, iNumDir), H(1, iNumDof);

          for (octave_idx_type i = 0; i < oIntegRule.iGetNumEvalPoints(); ++i) {
               const double alpha = oIntegRule.dGetWeight(i);

               for (octave_idx_type j = 0; j < iNumDir; ++j) {
                    rv.xelem(j) = oIntegRule.dGetPosition(i, j);
               }

               const double detJ = Jacobian(rv, J);

               AddMeshInfo(info, oIntegRule, detJ);

               ScalarInterpMatrix(rv, H, 0);

               for (octave_idx_type l = 0; l < iNumDof; ++l) {
                    for (octave_idx_type m = l; m < iNumDof; ++m) {
                         Ce.xelem(l, m) += H.xelem(l) * H.xelem(m) * alpha * rhocp * detJ;
                    }
               }
          }

          for (octave_idx_type i = 1; i < iNumDof; ++i) {
               for (octave_idx_type j = 0; j < i; ++j) {
                    Ce.xelem(i, j) = Ce.xelem(j, i);
               }
          }
     }
     
private:
     octave_idx_type iNumPreLoads;
     Matrix dTheta;
     NDArray epsilonRef;
};

class Iso8: public Element3D
{
public:
     Iso8(ElementTypes::TypeId eltype, octave_idx_type id, const Matrix& X, const Material* material, const int32NDArray& nodes, const StrainField& oRefStrain)
          :Element3D(eltype, id, X, material, nodes, oRefStrain) {
          FEM_ASSERT(nodes.numel() == 8);
     }

     virtual const IntegrationRule& GetIntegrationRule(FemMatrixType eMatType) const {
          static const octave_idx_type N = 2;
          static const double r[2][N] = {{0.577350269189626, -0.577350269189626}, {1., -1.}};
          static const double alpha[2][N] = {{1., 1.}, {1., 1.}};

          static array<IntegrationRule, 2> rgIntegRule;

          octave_idx_type iIntegRule;

          switch (eMatType) {
          case MAT_MASS_LUMPED:
               iIntegRule = 1;
               break;
          default:
               iIntegRule = 0;
          }

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

          return rgIntegRule[iIntegRule];
     }

protected:
     virtual double Jacobian(const ColumnVector& rv, Matrix& J) const {
          FEM_ASSERT(J.rows() == 3);
          FEM_ASSERT(J.columns() == 3);
          FEM_ASSERT(rv.numel() == 3);

          const double r = rv.xelem(0);
          const double s = rv.xelem(1);
          const double t = rv.xelem(2);

          J.xelem(0,0) = (-(X.xelem(0,1)*(s+1)*(t+1))/8.0)+(X.xelem(0,0)*(s+1)*(t+1))/8.0+(X.xelem(0,3)*(1-s)*(t+1))/8.0-(X.xelem(0,2)*(1-s)*(t+1))/8.0-(X.xelem(0,5)*(s+1)*(1-t))/8.0+(X.xelem(0,4)*(s+1)*(1-t))/8.0+(X.xelem(0,7)*(1-s)*(1-t))/8.0-(X.xelem(0,6)*(1-s)*(1-t))/8.0;
          J.xelem(1,0) = (-(X.xelem(0,3)*(r+1)*(t+1))/8.0)+(X.xelem(0,0)*(r+1)*(t+1))/8.0-(X.xelem(0,2)*(1-r)*(t+1))/8.0+(X.xelem(0,1)*(1-r)*(t+1))/8.0-(X.xelem(0,7)*(r+1)*(1-t))/8.0+(X.xelem(0,4)*(r+1)*(1-t))/8.0-(X.xelem(0,6)*(1-r)*(1-t))/8.0+(X.xelem(0,5)*(1-r)*(1-t))/8.0;
          J.xelem(2,0) = (-(X.xelem(0,4)*(r+1)*(s+1))/8.0)+(X.xelem(0,0)*(r+1)*(s+1))/8.0-(X.xelem(0,5)*(1-r)*(s+1))/8.0+(X.xelem(0,1)*(1-r)*(s+1))/8.0-(X.xelem(0,7)*(r+1)*(1-s))/8.0+(X.xelem(0,3)*(r+1)*(1-s))/8.0-(X.xelem(0,6)*(1-r)*(1-s))/8.0+(X.xelem(0,2)*(1-r)*(1-s))/8.0;
          J.xelem(0,1) = (-(X.xelem(1,1)*(s+1)*(t+1))/8.0)+(X.xelem(1,0)*(s+1)*(t+1))/8.0+(X.xelem(1,3)*(1-s)*(t+1))/8.0-(X.xelem(1,2)*(1-s)*(t+1))/8.0-(X.xelem(1,5)*(s+1)*(1-t))/8.0+(X.xelem(1,4)*(s+1)*(1-t))/8.0+(X.xelem(1,7)*(1-s)*(1-t))/8.0-(X.xelem(1,6)*(1-s)*(1-t))/8.0;
          J.xelem(1,1) = (-(X.xelem(1,3)*(r+1)*(t+1))/8.0)+(X.xelem(1,0)*(r+1)*(t+1))/8.0-(X.xelem(1,2)*(1-r)*(t+1))/8.0+(X.xelem(1,1)*(1-r)*(t+1))/8.0-(X.xelem(1,7)*(r+1)*(1-t))/8.0+(X.xelem(1,4)*(r+1)*(1-t))/8.0-(X.xelem(1,6)*(1-r)*(1-t))/8.0+(X.xelem(1,5)*(1-r)*(1-t))/8.0;
          J.xelem(2,1) = (-(X.xelem(1,4)*(r+1)*(s+1))/8.0)+(X.xelem(1,0)*(r+1)*(s+1))/8.0-(X.xelem(1,5)*(1-r)*(s+1))/8.0+(X.xelem(1,1)*(1-r)*(s+1))/8.0-(X.xelem(1,7)*(r+1)*(1-s))/8.0+(X.xelem(1,3)*(r+1)*(1-s))/8.0-(X.xelem(1,6)*(1-r)*(1-s))/8.0+(X.xelem(1,2)*(1-r)*(1-s))/8.0;
          J.xelem(0,2) = (-(X.xelem(2,1)*(s+1)*(t+1))/8.0)+(X.xelem(2,0)*(s+1)*(t+1))/8.0+(X.xelem(2,3)*(1-s)*(t+1))/8.0-(X.xelem(2,2)*(1-s)*(t+1))/8.0-(X.xelem(2,5)*(s+1)*(1-t))/8.0+(X.xelem(2,4)*(s+1)*(1-t))/8.0+(X.xelem(2,7)*(1-s)*(1-t))/8.0-(X.xelem(2,6)*(1-s)*(1-t))/8.0;
          J.xelem(1,2) = (-(X.xelem(2,3)*(r+1)*(t+1))/8.0)+(X.xelem(2,0)*(r+1)*(t+1))/8.0-(X.xelem(2,2)*(1-r)*(t+1))/8.0+(X.xelem(2,1)*(1-r)*(t+1))/8.0-(X.xelem(2,7)*(r+1)*(1-t))/8.0+(X.xelem(2,4)*(r+1)*(1-t))/8.0-(X.xelem(2,6)*(1-r)*(1-t))/8.0+(X.xelem(2,5)*(1-r)*(1-t))/8.0;
          J.xelem(2,2) = (-(X.xelem(2,4)*(r+1)*(s+1))/8.0)+(X.xelem(2,0)*(r+1)*(s+1))/8.0-(X.xelem(2,5)*(1-r)*(s+1))/8.0+(X.xelem(2,1)*(1-r)*(s+1))/8.0-(X.xelem(2,7)*(r+1)*(1-s))/8.0+(X.xelem(2,3)*(r+1)*(1-s))/8.0-(X.xelem(2,6)*(1-r)*(1-s))/8.0+(X.xelem(2,2)*(1-r)*(1-s))/8.0;

#ifdef DEBUG
          for (octave_idx_type i = 0; i < J.rows(); ++i) {
               for (octave_idx_type j = 0; j < J.columns(); ++j) {
                    FEM_ASSERT(std::isfinite(J(i, j)));
               }
          }
#endif
          return J.determinant();
     }

     virtual void DispInterpMatrix(const ColumnVector& rv, Matrix& H) const {
          FEM_ASSERT(rv.numel() == 3);
          FEM_ASSERT(H.rows() == 3);
          FEM_ASSERT(H.columns() == 24);

          const double r = rv.xelem(0);
          const double s = rv.xelem(1);
          const double t = rv.xelem(2);

          H.xelem(0,0) = ((r+1)*(s+1)*(t+1))/8.0;
          H.xelem(1,0) = 0;
          H.xelem(2,0) = 0;
          H.xelem(0,1) = 0;
          H.xelem(1,1) = ((r+1)*(s+1)*(t+1))/8.0;
          H.xelem(2,1) = 0;
          H.xelem(0,2) = 0;
          H.xelem(1,2) = 0;
          H.xelem(2,2) = ((r+1)*(s+1)*(t+1))/8.0;
          H.xelem(0,3) = ((1-r)*(s+1)*(t+1))/8.0;
          H.xelem(1,3) = 0;
          H.xelem(2,3) = 0;
          H.xelem(0,4) = 0;
          H.xelem(1,4) = ((1-r)*(s+1)*(t+1))/8.0;
          H.xelem(2,4) = 0;
          H.xelem(0,5) = 0;
          H.xelem(1,5) = 0;
          H.xelem(2,5) = ((1-r)*(s+1)*(t+1))/8.0;
          H.xelem(0,6) = ((1-r)*(1-s)*(t+1))/8.0;
          H.xelem(1,6) = 0;
          H.xelem(2,6) = 0;
          H.xelem(0,7) = 0;
          H.xelem(1,7) = ((1-r)*(1-s)*(t+1))/8.0;
          H.xelem(2,7) = 0;
          H.xelem(0,8) = 0;
          H.xelem(1,8) = 0;
          H.xelem(2,8) = ((1-r)*(1-s)*(t+1))/8.0;
          H.xelem(0,9) = ((r+1)*(1-s)*(t+1))/8.0;
          H.xelem(1,9) = 0;
          H.xelem(2,9) = 0;
          H.xelem(0,10) = 0;
          H.xelem(1,10) = ((r+1)*(1-s)*(t+1))/8.0;
          H.xelem(2,10) = 0;
          H.xelem(0,11) = 0;
          H.xelem(1,11) = 0;
          H.xelem(2,11) = ((r+1)*(1-s)*(t+1))/8.0;
          H.xelem(0,12) = ((r+1)*(s+1)*(1-t))/8.0;
          H.xelem(1,12) = 0;
          H.xelem(2,12) = 0;
          H.xelem(0,13) = 0;
          H.xelem(1,13) = ((r+1)*(s+1)*(1-t))/8.0;
          H.xelem(2,13) = 0;
          H.xelem(0,14) = 0;
          H.xelem(1,14) = 0;
          H.xelem(2,14) = ((r+1)*(s+1)*(1-t))/8.0;
          H.xelem(0,15) = ((1-r)*(s+1)*(1-t))/8.0;
          H.xelem(1,15) = 0;
          H.xelem(2,15) = 0;
          H.xelem(0,16) = 0;
          H.xelem(1,16) = ((1-r)*(s+1)*(1-t))/8.0;
          H.xelem(2,16) = 0;
          H.xelem(0,17) = 0;
          H.xelem(1,17) = 0;
          H.xelem(2,17) = ((1-r)*(s+1)*(1-t))/8.0;
          H.xelem(0,18) = ((1-r)*(1-s)*(1-t))/8.0;
          H.xelem(1,18) = 0;
          H.xelem(2,18) = 0;
          H.xelem(0,19) = 0;
          H.xelem(1,19) = ((1-r)*(1-s)*(1-t))/8.0;
          H.xelem(2,19) = 0;
          H.xelem(0,20) = 0;
          H.xelem(1,20) = 0;
          H.xelem(2,20) = ((1-r)*(1-s)*(1-t))/8.0;
          H.xelem(0,21) = ((r+1)*(1-s)*(1-t))/8.0;
          H.xelem(1,21) = 0;
          H.xelem(2,21) = 0;
          H.xelem(0,22) = 0;
          H.xelem(1,22) = ((r+1)*(1-s)*(1-t))/8.0;
          H.xelem(2,22) = 0;
          H.xelem(0,23) = 0;
          H.xelem(1,23) = 0;
          H.xelem(2,23) = ((r+1)*(1-s)*(1-t))/8.0;
     }

     virtual void StrainMatrix(const ColumnVector& rv, const Matrix& invJ, Matrix& B) const {
          FEM_ASSERT(rv.numel() == 3);
          FEM_ASSERT(invJ.rows() == 3);
          FEM_ASSERT(invJ.columns() == 3);
          FEM_ASSERT(B.rows() == 6);
          FEM_ASSERT(B.columns() == 24);

          const double r = rv.xelem(0);
          const double s = rv.xelem(1);
          const double t = rv.xelem(2);

          B.xelem(0,0) = (invJ.xelem(0,0)*(s+1)*(t+1))/8.0+(invJ.xelem(0,1)*(r+1)*(t+1))/8.0+(invJ.xelem(0,2)*(r+1)*(s+1))/8.0;
          B.xelem(1,0) = 0;
          B.xelem(2,0) = 0;
          B.xelem(3,0) = (invJ.xelem(1,0)*(s+1)*(t+1))/8.0+(invJ.xelem(1,1)*(r+1)*(t+1))/8.0+(invJ.xelem(1,2)*(r+1)*(s+1))/8.0;
          B.xelem(4,0) = 0;
          B.xelem(5,0) = (invJ.xelem(2,0)*(s+1)*(t+1))/8.0+(invJ.xelem(2,1)*(r+1)*(t+1))/8.0+(invJ.xelem(2,2)*(r+1)*(s+1))/8.0;
          B.xelem(0,1) = 0;
          B.xelem(1,1) = (invJ.xelem(1,0)*(s+1)*(t+1))/8.0+(invJ.xelem(1,1)*(r+1)*(t+1))/8.0+(invJ.xelem(1,2)*(r+1)*(s+1))/8.0;
          B.xelem(2,1) = 0;
          B.xelem(3,1) = (invJ.xelem(0,0)*(s+1)*(t+1))/8.0+(invJ.xelem(0,1)*(r+1)*(t+1))/8.0+(invJ.xelem(0,2)*(r+1)*(s+1))/8.0;
          B.xelem(4,1) = (invJ.xelem(2,0)*(s+1)*(t+1))/8.0+(invJ.xelem(2,1)*(r+1)*(t+1))/8.0+(invJ.xelem(2,2)*(r+1)*(s+1))/8.0;
          B.xelem(5,1) = 0;
          B.xelem(0,2) = 0;
          B.xelem(1,2) = 0;
          B.xelem(2,2) = (invJ.xelem(2,0)*(s+1)*(t+1))/8.0+(invJ.xelem(2,1)*(r+1)*(t+1))/8.0+(invJ.xelem(2,2)*(r+1)*(s+1))/8.0;
          B.xelem(3,2) = 0;
          B.xelem(4,2) = (invJ.xelem(1,0)*(s+1)*(t+1))/8.0+(invJ.xelem(1,1)*(r+1)*(t+1))/8.0+(invJ.xelem(1,2)*(r+1)*(s+1))/8.0;
          B.xelem(5,2) = (invJ.xelem(0,0)*(s+1)*(t+1))/8.0+(invJ.xelem(0,1)*(r+1)*(t+1))/8.0+(invJ.xelem(0,2)*(r+1)*(s+1))/8.0;
          B.xelem(0,3) = (-(invJ.xelem(0,0)*(s+1)*(t+1))/8.0)+(invJ.xelem(0,1)*(1-r)*(t+1))/8.0+(invJ.xelem(0,2)*(1-r)*(s+1))/8.0;
          B.xelem(1,3) = 0;
          B.xelem(2,3) = 0;
          B.xelem(3,3) = (-(invJ.xelem(1,0)*(s+1)*(t+1))/8.0)+(invJ.xelem(1,1)*(1-r)*(t+1))/8.0+(invJ.xelem(1,2)*(1-r)*(s+1))/8.0;
          B.xelem(4,3) = 0;
          B.xelem(5,3) = (-(invJ.xelem(2,0)*(s+1)*(t+1))/8.0)+(invJ.xelem(2,1)*(1-r)*(t+1))/8.0+(invJ.xelem(2,2)*(1-r)*(s+1))/8.0;
          B.xelem(0,4) = 0;
          B.xelem(1,4) = (-(invJ.xelem(1,0)*(s+1)*(t+1))/8.0)+(invJ.xelem(1,1)*(1-r)*(t+1))/8.0+(invJ.xelem(1,2)*(1-r)*(s+1))/8.0;
          B.xelem(2,4) = 0;
          B.xelem(3,4) = (-(invJ.xelem(0,0)*(s+1)*(t+1))/8.0)+(invJ.xelem(0,1)*(1-r)*(t+1))/8.0+(invJ.xelem(0,2)*(1-r)*(s+1))/8.0;
          B.xelem(4,4) = (-(invJ.xelem(2,0)*(s+1)*(t+1))/8.0)+(invJ.xelem(2,1)*(1-r)*(t+1))/8.0+(invJ.xelem(2,2)*(1-r)*(s+1))/8.0;
          B.xelem(5,4) = 0;
          B.xelem(0,5) = 0;
          B.xelem(1,5) = 0;
          B.xelem(2,5) = (-(invJ.xelem(2,0)*(s+1)*(t+1))/8.0)+(invJ.xelem(2,1)*(1-r)*(t+1))/8.0+(invJ.xelem(2,2)*(1-r)*(s+1))/8.0;
          B.xelem(3,5) = 0;
          B.xelem(4,5) = (-(invJ.xelem(1,0)*(s+1)*(t+1))/8.0)+(invJ.xelem(1,1)*(1-r)*(t+1))/8.0+(invJ.xelem(1,2)*(1-r)*(s+1))/8.0;
          B.xelem(5,5) = (-(invJ.xelem(0,0)*(s+1)*(t+1))/8.0)+(invJ.xelem(0,1)*(1-r)*(t+1))/8.0+(invJ.xelem(0,2)*(1-r)*(s+1))/8.0;
          B.xelem(0,6) = (-(invJ.xelem(0,0)*(1-s)*(t+1))/8.0)-(invJ.xelem(0,1)*(1-r)*(t+1))/8.0+(invJ.xelem(0,2)*(1-r)*(1-s))/8.0;
          B.xelem(1,6) = 0;
          B.xelem(2,6) = 0;
          B.xelem(3,6) = (-(invJ.xelem(1,0)*(1-s)*(t+1))/8.0)-(invJ.xelem(1,1)*(1-r)*(t+1))/8.0+(invJ.xelem(1,2)*(1-r)*(1-s))/8.0;
          B.xelem(4,6) = 0;
          B.xelem(5,6) = (-(invJ.xelem(2,0)*(1-s)*(t+1))/8.0)-(invJ.xelem(2,1)*(1-r)*(t+1))/8.0+(invJ.xelem(2,2)*(1-r)*(1-s))/8.0;
          B.xelem(0,7) = 0;
          B.xelem(1,7) = (-(invJ.xelem(1,0)*(1-s)*(t+1))/8.0)-(invJ.xelem(1,1)*(1-r)*(t+1))/8.0+(invJ.xelem(1,2)*(1-r)*(1-s))/8.0;
          B.xelem(2,7) = 0;
          B.xelem(3,7) = (-(invJ.xelem(0,0)*(1-s)*(t+1))/8.0)-(invJ.xelem(0,1)*(1-r)*(t+1))/8.0+(invJ.xelem(0,2)*(1-r)*(1-s))/8.0;
          B.xelem(4,7) = (-(invJ.xelem(2,0)*(1-s)*(t+1))/8.0)-(invJ.xelem(2,1)*(1-r)*(t+1))/8.0+(invJ.xelem(2,2)*(1-r)*(1-s))/8.0;
          B.xelem(5,7) = 0;
          B.xelem(0,8) = 0;
          B.xelem(1,8) = 0;
          B.xelem(2,8) = (-(invJ.xelem(2,0)*(1-s)*(t+1))/8.0)-(invJ.xelem(2,1)*(1-r)*(t+1))/8.0+(invJ.xelem(2,2)*(1-r)*(1-s))/8.0;
          B.xelem(3,8) = 0;
          B.xelem(4,8) = (-(invJ.xelem(1,0)*(1-s)*(t+1))/8.0)-(invJ.xelem(1,1)*(1-r)*(t+1))/8.0+(invJ.xelem(1,2)*(1-r)*(1-s))/8.0;
          B.xelem(5,8) = (-(invJ.xelem(0,0)*(1-s)*(t+1))/8.0)-(invJ.xelem(0,1)*(1-r)*(t+1))/8.0+(invJ.xelem(0,2)*(1-r)*(1-s))/8.0;
          B.xelem(0,9) = (invJ.xelem(0,0)*(1-s)*(t+1))/8.0-(invJ.xelem(0,1)*(r+1)*(t+1))/8.0+(invJ.xelem(0,2)*(r+1)*(1-s))/8.0;
          B.xelem(1,9) = 0;
          B.xelem(2,9) = 0;
          B.xelem(3,9) = (invJ.xelem(1,0)*(1-s)*(t+1))/8.0-(invJ.xelem(1,1)*(r+1)*(t+1))/8.0+(invJ.xelem(1,2)*(r+1)*(1-s))/8.0;
          B.xelem(4,9) = 0;
          B.xelem(5,9) = (invJ.xelem(2,0)*(1-s)*(t+1))/8.0-(invJ.xelem(2,1)*(r+1)*(t+1))/8.0+(invJ.xelem(2,2)*(r+1)*(1-s))/8.0;
          B.xelem(0,10) = 0;
          B.xelem(1,10) = (invJ.xelem(1,0)*(1-s)*(t+1))/8.0-(invJ.xelem(1,1)*(r+1)*(t+1))/8.0+(invJ.xelem(1,2)*(r+1)*(1-s))/8.0;
          B.xelem(2,10) = 0;
          B.xelem(3,10) = (invJ.xelem(0,0)*(1-s)*(t+1))/8.0-(invJ.xelem(0,1)*(r+1)*(t+1))/8.0+(invJ.xelem(0,2)*(r+1)*(1-s))/8.0;
          B.xelem(4,10) = (invJ.xelem(2,0)*(1-s)*(t+1))/8.0-(invJ.xelem(2,1)*(r+1)*(t+1))/8.0+(invJ.xelem(2,2)*(r+1)*(1-s))/8.0;
          B.xelem(5,10) = 0;
          B.xelem(0,11) = 0;
          B.xelem(1,11) = 0;
          B.xelem(2,11) = (invJ.xelem(2,0)*(1-s)*(t+1))/8.0-(invJ.xelem(2,1)*(r+1)*(t+1))/8.0+(invJ.xelem(2,2)*(r+1)*(1-s))/8.0;
          B.xelem(3,11) = 0;
          B.xelem(4,11) = (invJ.xelem(1,0)*(1-s)*(t+1))/8.0-(invJ.xelem(1,1)*(r+1)*(t+1))/8.0+(invJ.xelem(1,2)*(r+1)*(1-s))/8.0;
          B.xelem(5,11) = (invJ.xelem(0,0)*(1-s)*(t+1))/8.0-(invJ.xelem(0,1)*(r+1)*(t+1))/8.0+(invJ.xelem(0,2)*(r+1)*(1-s))/8.0;
          B.xelem(0,12) = (invJ.xelem(0,0)*(s+1)*(1-t))/8.0+(invJ.xelem(0,1)*(r+1)*(1-t))/8.0-(invJ.xelem(0,2)*(r+1)*(s+1))/8.0;
          B.xelem(1,12) = 0;
          B.xelem(2,12) = 0;
          B.xelem(3,12) = (invJ.xelem(1,0)*(s+1)*(1-t))/8.0+(invJ.xelem(1,1)*(r+1)*(1-t))/8.0-(invJ.xelem(1,2)*(r+1)*(s+1))/8.0;
          B.xelem(4,12) = 0;
          B.xelem(5,12) = (invJ.xelem(2,0)*(s+1)*(1-t))/8.0+(invJ.xelem(2,1)*(r+1)*(1-t))/8.0-(invJ.xelem(2,2)*(r+1)*(s+1))/8.0;
          B.xelem(0,13) = 0;
          B.xelem(1,13) = (invJ.xelem(1,0)*(s+1)*(1-t))/8.0+(invJ.xelem(1,1)*(r+1)*(1-t))/8.0-(invJ.xelem(1,2)*(r+1)*(s+1))/8.0;
          B.xelem(2,13) = 0;
          B.xelem(3,13) = (invJ.xelem(0,0)*(s+1)*(1-t))/8.0+(invJ.xelem(0,1)*(r+1)*(1-t))/8.0-(invJ.xelem(0,2)*(r+1)*(s+1))/8.0;
          B.xelem(4,13) = (invJ.xelem(2,0)*(s+1)*(1-t))/8.0+(invJ.xelem(2,1)*(r+1)*(1-t))/8.0-(invJ.xelem(2,2)*(r+1)*(s+1))/8.0;
          B.xelem(5,13) = 0;
          B.xelem(0,14) = 0;
          B.xelem(1,14) = 0;
          B.xelem(2,14) = (invJ.xelem(2,0)*(s+1)*(1-t))/8.0+(invJ.xelem(2,1)*(r+1)*(1-t))/8.0-(invJ.xelem(2,2)*(r+1)*(s+1))/8.0;
          B.xelem(3,14) = 0;
          B.xelem(4,14) = (invJ.xelem(1,0)*(s+1)*(1-t))/8.0+(invJ.xelem(1,1)*(r+1)*(1-t))/8.0-(invJ.xelem(1,2)*(r+1)*(s+1))/8.0;
          B.xelem(5,14) = (invJ.xelem(0,0)*(s+1)*(1-t))/8.0+(invJ.xelem(0,1)*(r+1)*(1-t))/8.0-(invJ.xelem(0,2)*(r+1)*(s+1))/8.0;
          B.xelem(0,15) = (-(invJ.xelem(0,0)*(s+1)*(1-t))/8.0)+(invJ.xelem(0,1)*(1-r)*(1-t))/8.0-(invJ.xelem(0,2)*(1-r)*(s+1))/8.0;
          B.xelem(1,15) = 0;
          B.xelem(2,15) = 0;
          B.xelem(3,15) = (-(invJ.xelem(1,0)*(s+1)*(1-t))/8.0)+(invJ.xelem(1,1)*(1-r)*(1-t))/8.0-(invJ.xelem(1,2)*(1-r)*(s+1))/8.0;
          B.xelem(4,15) = 0;
          B.xelem(5,15) = (-(invJ.xelem(2,0)*(s+1)*(1-t))/8.0)+(invJ.xelem(2,1)*(1-r)*(1-t))/8.0-(invJ.xelem(2,2)*(1-r)*(s+1))/8.0;
          B.xelem(0,16) = 0;
          B.xelem(1,16) = (-(invJ.xelem(1,0)*(s+1)*(1-t))/8.0)+(invJ.xelem(1,1)*(1-r)*(1-t))/8.0-(invJ.xelem(1,2)*(1-r)*(s+1))/8.0;
          B.xelem(2,16) = 0;
          B.xelem(3,16) = (-(invJ.xelem(0,0)*(s+1)*(1-t))/8.0)+(invJ.xelem(0,1)*(1-r)*(1-t))/8.0-(invJ.xelem(0,2)*(1-r)*(s+1))/8.0;
          B.xelem(4,16) = (-(invJ.xelem(2,0)*(s+1)*(1-t))/8.0)+(invJ.xelem(2,1)*(1-r)*(1-t))/8.0-(invJ.xelem(2,2)*(1-r)*(s+1))/8.0;
          B.xelem(5,16) = 0;
          B.xelem(0,17) = 0;
          B.xelem(1,17) = 0;
          B.xelem(2,17) = (-(invJ.xelem(2,0)*(s+1)*(1-t))/8.0)+(invJ.xelem(2,1)*(1-r)*(1-t))/8.0-(invJ.xelem(2,2)*(1-r)*(s+1))/8.0;
          B.xelem(3,17) = 0;
          B.xelem(4,17) = (-(invJ.xelem(1,0)*(s+1)*(1-t))/8.0)+(invJ.xelem(1,1)*(1-r)*(1-t))/8.0-(invJ.xelem(1,2)*(1-r)*(s+1))/8.0;
          B.xelem(5,17) = (-(invJ.xelem(0,0)*(s+1)*(1-t))/8.0)+(invJ.xelem(0,1)*(1-r)*(1-t))/8.0-(invJ.xelem(0,2)*(1-r)*(s+1))/8.0;
          B.xelem(0,18) = (-(invJ.xelem(0,0)*(1-s)*(1-t))/8.0)-(invJ.xelem(0,1)*(1-r)*(1-t))/8.0-(invJ.xelem(0,2)*(1-r)*(1-s))/8.0;
          B.xelem(1,18) = 0;
          B.xelem(2,18) = 0;
          B.xelem(3,18) = (-(invJ.xelem(1,0)*(1-s)*(1-t))/8.0)-(invJ.xelem(1,1)*(1-r)*(1-t))/8.0-(invJ.xelem(1,2)*(1-r)*(1-s))/8.0;
          B.xelem(4,18) = 0;
          B.xelem(5,18) = (-(invJ.xelem(2,0)*(1-s)*(1-t))/8.0)-(invJ.xelem(2,1)*(1-r)*(1-t))/8.0-(invJ.xelem(2,2)*(1-r)*(1-s))/8.0;
          B.xelem(0,19) = 0;
          B.xelem(1,19) = (-(invJ.xelem(1,0)*(1-s)*(1-t))/8.0)-(invJ.xelem(1,1)*(1-r)*(1-t))/8.0-(invJ.xelem(1,2)*(1-r)*(1-s))/8.0;
          B.xelem(2,19) = 0;
          B.xelem(3,19) = (-(invJ.xelem(0,0)*(1-s)*(1-t))/8.0)-(invJ.xelem(0,1)*(1-r)*(1-t))/8.0-(invJ.xelem(0,2)*(1-r)*(1-s))/8.0;
          B.xelem(4,19) = (-(invJ.xelem(2,0)*(1-s)*(1-t))/8.0)-(invJ.xelem(2,1)*(1-r)*(1-t))/8.0-(invJ.xelem(2,2)*(1-r)*(1-s))/8.0;
          B.xelem(5,19) = 0;
          B.xelem(0,20) = 0;
          B.xelem(1,20) = 0;
          B.xelem(2,20) = (-(invJ.xelem(2,0)*(1-s)*(1-t))/8.0)-(invJ.xelem(2,1)*(1-r)*(1-t))/8.0-(invJ.xelem(2,2)*(1-r)*(1-s))/8.0;
          B.xelem(3,20) = 0;
          B.xelem(4,20) = (-(invJ.xelem(1,0)*(1-s)*(1-t))/8.0)-(invJ.xelem(1,1)*(1-r)*(1-t))/8.0-(invJ.xelem(1,2)*(1-r)*(1-s))/8.0;
          B.xelem(5,20) = (-(invJ.xelem(0,0)*(1-s)*(1-t))/8.0)-(invJ.xelem(0,1)*(1-r)*(1-t))/8.0-(invJ.xelem(0,2)*(1-r)*(1-s))/8.0;
          B.xelem(0,21) = (invJ.xelem(0,0)*(1-s)*(1-t))/8.0-(invJ.xelem(0,1)*(r+1)*(1-t))/8.0-(invJ.xelem(0,2)*(r+1)*(1-s))/8.0;
          B.xelem(1,21) = 0;
          B.xelem(2,21) = 0;
          B.xelem(3,21) = (invJ.xelem(1,0)*(1-s)*(1-t))/8.0-(invJ.xelem(1,1)*(r+1)*(1-t))/8.0-(invJ.xelem(1,2)*(r+1)*(1-s))/8.0;
          B.xelem(4,21) = 0;
          B.xelem(5,21) = (invJ.xelem(2,0)*(1-s)*(1-t))/8.0-(invJ.xelem(2,1)*(r+1)*(1-t))/8.0-(invJ.xelem(2,2)*(r+1)*(1-s))/8.0;
          B.xelem(0,22) = 0;
          B.xelem(1,22) = (invJ.xelem(1,0)*(1-s)*(1-t))/8.0-(invJ.xelem(1,1)*(r+1)*(1-t))/8.0-(invJ.xelem(1,2)*(r+1)*(1-s))/8.0;
          B.xelem(2,22) = 0;
          B.xelem(3,22) = (invJ.xelem(0,0)*(1-s)*(1-t))/8.0-(invJ.xelem(0,1)*(r+1)*(1-t))/8.0-(invJ.xelem(0,2)*(r+1)*(1-s))/8.0;
          B.xelem(4,22) = (invJ.xelem(2,0)*(1-s)*(1-t))/8.0-(invJ.xelem(2,1)*(r+1)*(1-t))/8.0-(invJ.xelem(2,2)*(r+1)*(1-s))/8.0;
          B.xelem(5,22) = 0;
          B.xelem(0,23) = 0;
          B.xelem(1,23) = 0;
          B.xelem(2,23) = (invJ.xelem(2,0)*(1-s)*(1-t))/8.0-(invJ.xelem(2,1)*(r+1)*(1-t))/8.0-(invJ.xelem(2,2)*(r+1)*(1-s))/8.0;
          B.xelem(3,23) = 0;
          B.xelem(4,23) = (invJ.xelem(1,0)*(1-s)*(1-t))/8.0-(invJ.xelem(1,1)*(r+1)*(1-t))/8.0-(invJ.xelem(1,2)*(r+1)*(1-s))/8.0;
          B.xelem(5,23) = (invJ.xelem(0,0)*(1-s)*(1-t))/8.0-(invJ.xelem(0,1)*(r+1)*(1-t))/8.0-(invJ.xelem(0,2)*(r+1)*(1-s))/8.0;

#ifdef DEBUG
          for (octave_idx_type i = 0; i < B.rows(); ++i) {
               for (octave_idx_type j = 0; j < B.columns(); ++j) {
                    FEM_ASSERT(std::isfinite(B.xelem(i, j)));
               }
          }
#endif
     }

     virtual void TemperatureGradientMatrix(const ColumnVector& rv, const Matrix& invJ, Matrix& Bt) const {
          FEM_ASSERT(rv.numel() == 3);
          FEM_ASSERT(invJ.rows() == 3);
          FEM_ASSERT(invJ.columns() == 3);
          FEM_ASSERT(Bt.rows() == 3);
          FEM_ASSERT(Bt.columns() == 8);

          const double r = rv.xelem(0);
          const double s = rv.xelem(1);
          const double t = rv.xelem(2);
          
          // temperature gradient matrix Bt
          Bt.xelem(0,0) = (invJ.xelem(0,0)*(s+1)*(t+1))/8.0E+0+(invJ.xelem(0,1)*(r+1)*(t+1))/8.0E+0+(invJ.xelem(0,2)*(r+1)*(s+1))/8.0E+0;
          Bt.xelem(1,0) = (invJ.xelem(1,0)*(s+1)*(t+1))/8.0E+0+(invJ.xelem(1,1)*(r+1)*(t+1))/8.0E+0+(invJ.xelem(1,2)*(r+1)*(s+1))/8.0E+0;
          Bt.xelem(2,0) = (invJ.xelem(2,0)*(s+1)*(t+1))/8.0E+0+(invJ.xelem(2,1)*(r+1)*(t+1))/8.0E+0+(invJ.xelem(2,2)*(r+1)*(s+1))/8.0E+0;
          Bt.xelem(0,1) = (-(invJ.xelem(0,0)*(s+1)*(t+1))/8.0E+0)+(invJ.xelem(0,1)*(1-r)*(t+1))/8.0E+0+(invJ.xelem(0,2)*(1-r)*(s+1))/8.0E+0;
          Bt.xelem(1,1) = (-(invJ.xelem(1,0)*(s+1)*(t+1))/8.0E+0)+(invJ.xelem(1,1)*(1-r)*(t+1))/8.0E+0+(invJ.xelem(1,2)*(1-r)*(s+1))/8.0E+0;
          Bt.xelem(2,1) = (-(invJ.xelem(2,0)*(s+1)*(t+1))/8.0E+0)+(invJ.xelem(2,1)*(1-r)*(t+1))/8.0E+0+(invJ.xelem(2,2)*(1-r)*(s+1))/8.0E+0;
          Bt.xelem(0,2) = (-(invJ.xelem(0,0)*(1-s)*(t+1))/8.0E+0)-(invJ.xelem(0,1)*(1-r)*(t+1))/8.0E+0+(invJ.xelem(0,2)*(1-r)*(1-s))/8.0E+0;
          Bt.xelem(1,2) = (-(invJ.xelem(1,0)*(1-s)*(t+1))/8.0E+0)-(invJ.xelem(1,1)*(1-r)*(t+1))/8.0E+0+(invJ.xelem(1,2)*(1-r)*(1-s))/8.0E+0;
          Bt.xelem(2,2) = (-(invJ.xelem(2,0)*(1-s)*(t+1))/8.0E+0)-(invJ.xelem(2,1)*(1-r)*(t+1))/8.0E+0+(invJ.xelem(2,2)*(1-r)*(1-s))/8.0E+0;
          Bt.xelem(0,3) = (invJ.xelem(0,0)*(1-s)*(t+1))/8.0E+0-(invJ.xelem(0,1)*(r+1)*(t+1))/8.0E+0+(invJ.xelem(0,2)*(r+1)*(1-s))/8.0E+0;
          Bt.xelem(1,3) = (invJ.xelem(1,0)*(1-s)*(t+1))/8.0E+0-(invJ.xelem(1,1)*(r+1)*(t+1))/8.0E+0+(invJ.xelem(1,2)*(r+1)*(1-s))/8.0E+0;
          Bt.xelem(2,3) = (invJ.xelem(2,0)*(1-s)*(t+1))/8.0E+0-(invJ.xelem(2,1)*(r+1)*(t+1))/8.0E+0+(invJ.xelem(2,2)*(r+1)*(1-s))/8.0E+0;
          Bt.xelem(0,4) = (invJ.xelem(0,0)*(s+1)*(1-t))/8.0E+0+(invJ.xelem(0,1)*(r+1)*(1-t))/8.0E+0-(invJ.xelem(0,2)*(r+1)*(s+1))/8.0E+0;
          Bt.xelem(1,4) = (invJ.xelem(1,0)*(s+1)*(1-t))/8.0E+0+(invJ.xelem(1,1)*(r+1)*(1-t))/8.0E+0-(invJ.xelem(1,2)*(r+1)*(s+1))/8.0E+0;
          Bt.xelem(2,4) = (invJ.xelem(2,0)*(s+1)*(1-t))/8.0E+0+(invJ.xelem(2,1)*(r+1)*(1-t))/8.0E+0-(invJ.xelem(2,2)*(r+1)*(s+1))/8.0E+0;
          Bt.xelem(0,5) = (-(invJ.xelem(0,0)*(s+1)*(1-t))/8.0E+0)+(invJ.xelem(0,1)*(1-r)*(1-t))/8.0E+0-(invJ.xelem(0,2)*(1-r)*(s+1))/8.0E+0;
          Bt.xelem(1,5) = (-(invJ.xelem(1,0)*(s+1)*(1-t))/8.0E+0)+(invJ.xelem(1,1)*(1-r)*(1-t))/8.0E+0-(invJ.xelem(1,2)*(1-r)*(s+1))/8.0E+0;
          Bt.xelem(2,5) = (-(invJ.xelem(2,0)*(s+1)*(1-t))/8.0E+0)+(invJ.xelem(2,1)*(1-r)*(1-t))/8.0E+0-(invJ.xelem(2,2)*(1-r)*(s+1))/8.0E+0;
          Bt.xelem(0,6) = (-(invJ.xelem(0,0)*(1-s)*(1-t))/8.0E+0)-(invJ.xelem(0,1)*(1-r)*(1-t))/8.0E+0-(invJ.xelem(0,2)*(1-r)*(1-s))/8.0E+0;
          Bt.xelem(1,6) = (-(invJ.xelem(1,0)*(1-s)*(1-t))/8.0E+0)-(invJ.xelem(1,1)*(1-r)*(1-t))/8.0E+0-(invJ.xelem(1,2)*(1-r)*(1-s))/8.0E+0;
          Bt.xelem(2,6) = (-(invJ.xelem(2,0)*(1-s)*(1-t))/8.0E+0)-(invJ.xelem(2,1)*(1-r)*(1-t))/8.0E+0-(invJ.xelem(2,2)*(1-r)*(1-s))/8.0E+0;
          Bt.xelem(0,7) = (invJ.xelem(0,0)*(1-s)*(1-t))/8.0E+0-(invJ.xelem(0,1)*(r+1)*(1-t))/8.0E+0-(invJ.xelem(0,2)*(r+1)*(1-s))/8.0E+0;
          Bt.xelem(1,7) = (invJ.xelem(1,0)*(1-s)*(1-t))/8.0E+0-(invJ.xelem(1,1)*(r+1)*(1-t))/8.0E+0-(invJ.xelem(1,2)*(r+1)*(1-s))/8.0E+0;
          Bt.xelem(2,7) = (invJ.xelem(2,0)*(1-s)*(1-t))/8.0E+0-(invJ.xelem(2,1)*(r+1)*(1-t))/8.0E+0-(invJ.xelem(2,2)*(r+1)*(1-s))/8.0E+0;
     }
     
     virtual Matrix InterpGaussToNodal(FemMatrixType eMatType, const Matrix& taug) const {
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

private:
     void ScalarInterpMatrix(const ColumnVector& rv, Matrix& Hs, octave_idx_type irow) const {
          FEM_ASSERT(rv.numel() == 3);
          FEM_ASSERT(Hs.columns() == 8);
          FEM_ASSERT(irow >= 0);
          FEM_ASSERT(irow < Hs.rows());

          const double r = rv.xelem(0);
          const double s = rv.xelem(1);
          const double t = rv.xelem(2);

          Hs.xelem(irow, 0) = ((r+1)*(s+1)*(t+1))/8.0;
          Hs.xelem(irow, 1) = ((1-r)*(s+1)*(t+1))/8.0;
          Hs.xelem(irow, 2) = ((1-r)*(1-s)*(t+1))/8.0;
          Hs.xelem(irow, 3) = ((r+1)*(1-s)*(t+1))/8.0;
          Hs.xelem(irow, 4) = ((r+1)*(s+1)*(1-t))/8.0;
          Hs.xelem(irow, 5) = ((1-r)*(s+1)*(1-t))/8.0;
          Hs.xelem(irow, 6) = ((1-r)*(1-s)*(1-t))/8.0;
          Hs.xelem(irow, 7) = ((r+1)*(1-s)*(1-t))/8.0;
     }
};


class Iso20: public Element3D
{
public:
     Iso20(ElementTypes::TypeId eltype, octave_idx_type id, const Matrix& X, const Material* material, const int32NDArray& nodes, const StrainField& oRefStrain)
          :Element3D(eltype, id, X, material, nodes, oRefStrain) {
          FEM_ASSERT(nodes.numel() == 20);
     }

     virtual const IntegrationRule& GetIntegrationRule(FemMatrixType eMatType) const {
          constexpr octave_idx_type N = 3;
          static const double r[2][N] = {{0.774596669241483, 0., -0.774596669241483}, {1., 0., -1.}};
          static const double alpha[2][N] = {{0.555555555555556, 0.888888888888889, 0.555555555555556}, {2./3., 2./3., 2./3.}};

          static array<IntegrationRule, 2> rgIntegRule;

          octave_idx_type iIntegRule;

          switch (eMatType) {
          case MAT_MASS_LUMPED:
               iIntegRule = 1;
               break;
          default:
               iIntegRule = 0;
          }

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

          return rgIntegRule[iIntegRule];
     }

protected:
     virtual double Jacobian(const ColumnVector& rv, Matrix& J) const {
          FEM_ASSERT(J.rows() == 3);
          FEM_ASSERT(J.columns() == 3);
          FEM_ASSERT(rv.numel() == 3);

          const double r = rv.xelem(0);
          const double s = rv.xelem(1);
          const double t = rv.xelem(2);
          const double r2 = r * r;
          const double s2 = s * s;
          const double t2 = t * t;

          J.xelem(0,0) = (-(X.xelem(0,17)*(s+1)*(1-t2))/4.0E+0)+(X.xelem(0,16)*(s+1)*(1-t2))/4.0E+0+(X.xelem(0,19)*(1-s)*(1-t2))/4.0E+0-(X.xelem(0,18)*(1-s)*(1-t2))/4.0E+0+X.xelem(0,0)*(((s+1)*(t+1))/8.0E+0-(((s+1)*(1-t2))/4.0E+0+((1-s2)*(t+1))/4.0E+0-(r*(s+1)*(t+1))/2.0E+0)/2.0E+0)+X.xelem(0,4)*(((s+1)*(1-t))/8.0E+0-(((s+1)*(1-t2))/4.0E+0+((1-s2)*(1-t))/4.0E+0-(r*(s+1)*(1-t))/2.0E+0)/2.0E+0)+X.xelem(0,1)*((-((-((s+1)*(1-t2))/4.0E+0)-((1-s2)*(t+1))/4.0E+0-(r*(s+1)*(t+1))/2.0E+0)/2.0E+0)-((s+1)*(t+1))/8.0E+0)+X.xelem(0,5)*((-((-((s+1)*(1-t2))/4.0E+0)-((1-s2)*(1-t))/4.0E+0-(r*(s+1)*(1-t))/2.0E+0)/2.0E+0)-((s+1)*(1-t))/8.0E+0)+X.xelem(0,3)*(((1-s)*(t+1))/8.0E+0-(((1-s)*(1-t2))/4.0E+0+((1-s2)*(t+1))/4.0E+0-(r*(1-s)*(t+1))/2.0E+0)/2.0E+0)+X.xelem(0,7)*(((1-s)*(1-t))/8.0E+0-(((1-s)*(1-t2))/4.0E+0+((1-s2)*(1-t))/4.0E+0-(r*(1-s)*(1-t))/2.0E+0)/2.0E+0)+X.xelem(0,2)*((-((-((1-s)*(1-t2))/4.0E+0)-((1-s2)*(t+1))/4.0E+0-(r*(1-s)*(t+1))/2.0E+0)/2.0E+0)-((1-s)*(t+1))/8.0E+0)+X.xelem(0,6)*((-((-((1-s)*(1-t2))/4.0E+0)-((1-s2)*(1-t))/4.0E+0-(r*(1-s)*(1-t))/2.0E+0)/2.0E+0)-((1-s)*(1-t))/8.0E+0)+(X.xelem(0,11)*(1-s2)*(t+1))/4.0E+0-(X.xelem(0,9)*(1-s2)*(t+1))/4.0E+0-(X.xelem(0,8)*r*(s+1)*(t+1))/2.0E+0-(X.xelem(0,10)*r*(1-s)*(t+1))/2.0E+0+(X.xelem(0,15)*(1-s2)*(1-t))/4.0E+0-(X.xelem(0,13)*(1-s2)*(1-t))/4.0E+0-(X.xelem(0,12)*r*(s+1)*(1-t))/2.0E+0-(X.xelem(0,14)*r*(1-s)*(1-t))/2.0E+0;
          J.xelem(1,0) = (-(X.xelem(0,19)*(r+1)*(1-t2))/4.0E+0)+(X.xelem(0,16)*(r+1)*(1-t2))/4.0E+0-(X.xelem(0,18)*(1-r)*(1-t2))/4.0E+0+(X.xelem(0,17)*(1-r)*(1-t2))/4.0E+0+X.xelem(0,0)*(((r+1)*(t+1))/8.0E+0-(((r+1)*(1-t2))/4.0E+0-((r+1)*s*(t+1))/2.0E+0+((1-r2)*(t+1))/4.0E+0)/2.0E+0)+X.xelem(0,4)*(((r+1)*(1-t))/8.0E+0-(((r+1)*(1-t2))/4.0E+0-((r+1)*s*(1-t))/2.0E+0+((1-r2)*(1-t))/4.0E+0)/2.0E+0)+X.xelem(0,3)*((-((-((r+1)*(1-t2))/4.0E+0)-((r+1)*s*(t+1))/2.0E+0-((1-r2)*(t+1))/4.0E+0)/2.0E+0)-((r+1)*(t+1))/8.0E+0)+X.xelem(0,7)*((-((-((r+1)*(1-t2))/4.0E+0)-((r+1)*s*(1-t))/2.0E+0-((1-r2)*(1-t))/4.0E+0)/2.0E+0)-((r+1)*(1-t))/8.0E+0)+X.xelem(0,1)*(((1-r)*(t+1))/8.0E+0-(((1-r)*(1-t2))/4.0E+0-((1-r)*s*(t+1))/2.0E+0+((1-r2)*(t+1))/4.0E+0)/2.0E+0)+X.xelem(0,5)*(((1-r)*(1-t))/8.0E+0-(((1-r)*(1-t2))/4.0E+0-((1-r)*s*(1-t))/2.0E+0+((1-r2)*(1-t))/4.0E+0)/2.0E+0)+X.xelem(0,2)*((-((-((1-r)*(1-t2))/4.0E+0)-((1-r)*s*(t+1))/2.0E+0-((1-r2)*(t+1))/4.0E+0)/2.0E+0)-((1-r)*(t+1))/8.0E+0)+X.xelem(0,6)*((-((-((1-r)*(1-t2))/4.0E+0)-((1-r)*s*(1-t))/2.0E+0-((1-r2)*(1-t))/4.0E+0)/2.0E+0)-((1-r)*(1-t))/8.0E+0)-(X.xelem(0,11)*(r+1)*s*(t+1))/2.0E+0-(X.xelem(0,9)*(1-r)*s*(t+1))/2.0E+0-(X.xelem(0,10)*(1-r2)*(t+1))/4.0E+0+(X.xelem(0,8)*(1-r2)*(t+1))/4.0E+0-(X.xelem(0,15)*(r+1)*s*(1-t))/2.0E+0-(X.xelem(0,13)*(1-r)*s*(1-t))/2.0E+0-(X.xelem(0,14)*(1-r2)*(1-t))/4.0E+0+(X.xelem(0,12)*(1-r2)*(1-t))/4.0E+0;
          J.xelem(2,0) = X.xelem(0,0)*(((r+1)*(s+1))/8.0E+0-((-((r+1)*(s+1)*t)/2.0E+0)+((r+1)*(1-s2))/4.0E+0+((1-r2)*(s+1))/4.0E+0)/2.0E+0)+X.xelem(0,4)*((-((-((r+1)*(s+1)*t)/2.0E+0)-((r+1)*(1-s2))/4.0E+0-((1-r2)*(s+1))/4.0E+0)/2.0E+0)-((r+1)*(s+1))/8.0E+0)+X.xelem(0,1)*(((1-r)*(s+1))/8.0E+0-((-((1-r)*(s+1)*t)/2.0E+0)+((1-r)*(1-s2))/4.0E+0+((1-r2)*(s+1))/4.0E+0)/2.0E+0)+X.xelem(0,5)*((-((-((1-r)*(s+1)*t)/2.0E+0)-((1-r)*(1-s2))/4.0E+0-((1-r2)*(s+1))/4.0E+0)/2.0E+0)-((1-r)*(s+1))/8.0E+0)+X.xelem(0,3)*(((r+1)*(1-s))/8.0E+0-((-((r+1)*(1-s)*t)/2.0E+0)+((r+1)*(1-s2))/4.0E+0+((1-r2)*(1-s))/4.0E+0)/2.0E+0)+X.xelem(0,7)*((-((-((r+1)*(1-s)*t)/2.0E+0)-((r+1)*(1-s2))/4.0E+0-((1-r2)*(1-s))/4.0E+0)/2.0E+0)-((r+1)*(1-s))/8.0E+0)+X.xelem(0,2)*(((1-r)*(1-s))/8.0E+0-((-((1-r)*(1-s)*t)/2.0E+0)+((1-r)*(1-s2))/4.0E+0+((1-r2)*(1-s))/4.0E+0)/2.0E+0)+X.xelem(0,6)*((-((-((1-r)*(1-s)*t)/2.0E+0)-((1-r)*(1-s2))/4.0E+0-((1-r2)*(1-s))/4.0E+0)/2.0E+0)-((1-r)*(1-s))/8.0E+0)-(X.xelem(0,16)*(r+1)*(s+1)*t)/2.0E+0-(X.xelem(0,17)*(1-r)*(s+1)*t)/2.0E+0-(X.xelem(0,19)*(r+1)*(1-s)*t)/2.0E+0-(X.xelem(0,18)*(1-r)*(1-s)*t)/2.0E+0-(X.xelem(0,15)*(r+1)*(1-s2))/4.0E+0+(X.xelem(0,11)*(r+1)*(1-s2))/4.0E+0-(X.xelem(0,13)*(1-r)*(1-s2))/4.0E+0+(X.xelem(0,9)*(1-r)*(1-s2))/4.0E+0-(X.xelem(0,12)*(1-r2)*(s+1))/4.0E+0+(X.xelem(0,8)*(1-r2)*(s+1))/4.0E+0-(X.xelem(0,14)*(1-r2)*(1-s))/4.0E+0+(X.xelem(0,10)*(1-r2)*(1-s))/4.0E+0;
          J.xelem(0,1) = (-(X.xelem(1,17)*(s+1)*(1-t2))/4.0E+0)+(X.xelem(1,16)*(s+1)*(1-t2))/4.0E+0+(X.xelem(1,19)*(1-s)*(1-t2))/4.0E+0-(X.xelem(1,18)*(1-s)*(1-t2))/4.0E+0+X.xelem(1,0)*(((s+1)*(t+1))/8.0E+0-(((s+1)*(1-t2))/4.0E+0+((1-s2)*(t+1))/4.0E+0-(r*(s+1)*(t+1))/2.0E+0)/2.0E+0)+X.xelem(1,4)*(((s+1)*(1-t))/8.0E+0-(((s+1)*(1-t2))/4.0E+0+((1-s2)*(1-t))/4.0E+0-(r*(s+1)*(1-t))/2.0E+0)/2.0E+0)+X.xelem(1,1)*((-((-((s+1)*(1-t2))/4.0E+0)-((1-s2)*(t+1))/4.0E+0-(r*(s+1)*(t+1))/2.0E+0)/2.0E+0)-((s+1)*(t+1))/8.0E+0)+X.xelem(1,5)*((-((-((s+1)*(1-t2))/4.0E+0)-((1-s2)*(1-t))/4.0E+0-(r*(s+1)*(1-t))/2.0E+0)/2.0E+0)-((s+1)*(1-t))/8.0E+0)+X.xelem(1,3)*(((1-s)*(t+1))/8.0E+0-(((1-s)*(1-t2))/4.0E+0+((1-s2)*(t+1))/4.0E+0-(r*(1-s)*(t+1))/2.0E+0)/2.0E+0)+X.xelem(1,7)*(((1-s)*(1-t))/8.0E+0-(((1-s)*(1-t2))/4.0E+0+((1-s2)*(1-t))/4.0E+0-(r*(1-s)*(1-t))/2.0E+0)/2.0E+0)+X.xelem(1,2)*((-((-((1-s)*(1-t2))/4.0E+0)-((1-s2)*(t+1))/4.0E+0-(r*(1-s)*(t+1))/2.0E+0)/2.0E+0)-((1-s)*(t+1))/8.0E+0)+X.xelem(1,6)*((-((-((1-s)*(1-t2))/4.0E+0)-((1-s2)*(1-t))/4.0E+0-(r*(1-s)*(1-t))/2.0E+0)/2.0E+0)-((1-s)*(1-t))/8.0E+0)+(X.xelem(1,11)*(1-s2)*(t+1))/4.0E+0-(X.xelem(1,9)*(1-s2)*(t+1))/4.0E+0-(X.xelem(1,8)*r*(s+1)*(t+1))/2.0E+0-(X.xelem(1,10)*r*(1-s)*(t+1))/2.0E+0+(X.xelem(1,15)*(1-s2)*(1-t))/4.0E+0-(X.xelem(1,13)*(1-s2)*(1-t))/4.0E+0-(X.xelem(1,12)*r*(s+1)*(1-t))/2.0E+0-(X.xelem(1,14)*r*(1-s)*(1-t))/2.0E+0;
          J.xelem(1,1) = (-(X.xelem(1,19)*(r+1)*(1-t2))/4.0E+0)+(X.xelem(1,16)*(r+1)*(1-t2))/4.0E+0-(X.xelem(1,18)*(1-r)*(1-t2))/4.0E+0+(X.xelem(1,17)*(1-r)*(1-t2))/4.0E+0+X.xelem(1,0)*(((r+1)*(t+1))/8.0E+0-(((r+1)*(1-t2))/4.0E+0-((r+1)*s*(t+1))/2.0E+0+((1-r2)*(t+1))/4.0E+0)/2.0E+0)+X.xelem(1,4)*(((r+1)*(1-t))/8.0E+0-(((r+1)*(1-t2))/4.0E+0-((r+1)*s*(1-t))/2.0E+0+((1-r2)*(1-t))/4.0E+0)/2.0E+0)+X.xelem(1,3)*((-((-((r+1)*(1-t2))/4.0E+0)-((r+1)*s*(t+1))/2.0E+0-((1-r2)*(t+1))/4.0E+0)/2.0E+0)-((r+1)*(t+1))/8.0E+0)+X.xelem(1,7)*((-((-((r+1)*(1-t2))/4.0E+0)-((r+1)*s*(1-t))/2.0E+0-((1-r2)*(1-t))/4.0E+0)/2.0E+0)-((r+1)*(1-t))/8.0E+0)+X.xelem(1,1)*(((1-r)*(t+1))/8.0E+0-(((1-r)*(1-t2))/4.0E+0-((1-r)*s*(t+1))/2.0E+0+((1-r2)*(t+1))/4.0E+0)/2.0E+0)+X.xelem(1,5)*(((1-r)*(1-t))/8.0E+0-(((1-r)*(1-t2))/4.0E+0-((1-r)*s*(1-t))/2.0E+0+((1-r2)*(1-t))/4.0E+0)/2.0E+0)+X.xelem(1,2)*((-((-((1-r)*(1-t2))/4.0E+0)-((1-r)*s*(t+1))/2.0E+0-((1-r2)*(t+1))/4.0E+0)/2.0E+0)-((1-r)*(t+1))/8.0E+0)+X.xelem(1,6)*((-((-((1-r)*(1-t2))/4.0E+0)-((1-r)*s*(1-t))/2.0E+0-((1-r2)*(1-t))/4.0E+0)/2.0E+0)-((1-r)*(1-t))/8.0E+0)-(X.xelem(1,11)*(r+1)*s*(t+1))/2.0E+0-(X.xelem(1,9)*(1-r)*s*(t+1))/2.0E+0-(X.xelem(1,10)*(1-r2)*(t+1))/4.0E+0+(X.xelem(1,8)*(1-r2)*(t+1))/4.0E+0-(X.xelem(1,15)*(r+1)*s*(1-t))/2.0E+0-(X.xelem(1,13)*(1-r)*s*(1-t))/2.0E+0-(X.xelem(1,14)*(1-r2)*(1-t))/4.0E+0+(X.xelem(1,12)*(1-r2)*(1-t))/4.0E+0;
          J.xelem(2,1) = X.xelem(1,0)*(((r+1)*(s+1))/8.0E+0-((-((r+1)*(s+1)*t)/2.0E+0)+((r+1)*(1-s2))/4.0E+0+((1-r2)*(s+1))/4.0E+0)/2.0E+0)+X.xelem(1,4)*((-((-((r+1)*(s+1)*t)/2.0E+0)-((r+1)*(1-s2))/4.0E+0-((1-r2)*(s+1))/4.0E+0)/2.0E+0)-((r+1)*(s+1))/8.0E+0)+X.xelem(1,1)*(((1-r)*(s+1))/8.0E+0-((-((1-r)*(s+1)*t)/2.0E+0)+((1-r)*(1-s2))/4.0E+0+((1-r2)*(s+1))/4.0E+0)/2.0E+0)+X.xelem(1,5)*((-((-((1-r)*(s+1)*t)/2.0E+0)-((1-r)*(1-s2))/4.0E+0-((1-r2)*(s+1))/4.0E+0)/2.0E+0)-((1-r)*(s+1))/8.0E+0)+X.xelem(1,3)*(((r+1)*(1-s))/8.0E+0-((-((r+1)*(1-s)*t)/2.0E+0)+((r+1)*(1-s2))/4.0E+0+((1-r2)*(1-s))/4.0E+0)/2.0E+0)+X.xelem(1,7)*((-((-((r+1)*(1-s)*t)/2.0E+0)-((r+1)*(1-s2))/4.0E+0-((1-r2)*(1-s))/4.0E+0)/2.0E+0)-((r+1)*(1-s))/8.0E+0)+X.xelem(1,2)*(((1-r)*(1-s))/8.0E+0-((-((1-r)*(1-s)*t)/2.0E+0)+((1-r)*(1-s2))/4.0E+0+((1-r2)*(1-s))/4.0E+0)/2.0E+0)+X.xelem(1,6)*((-((-((1-r)*(1-s)*t)/2.0E+0)-((1-r)*(1-s2))/4.0E+0-((1-r2)*(1-s))/4.0E+0)/2.0E+0)-((1-r)*(1-s))/8.0E+0)-(X.xelem(1,16)*(r+1)*(s+1)*t)/2.0E+0-(X.xelem(1,17)*(1-r)*(s+1)*t)/2.0E+0-(X.xelem(1,19)*(r+1)*(1-s)*t)/2.0E+0-(X.xelem(1,18)*(1-r)*(1-s)*t)/2.0E+0-(X.xelem(1,15)*(r+1)*(1-s2))/4.0E+0+(X.xelem(1,11)*(r+1)*(1-s2))/4.0E+0-(X.xelem(1,13)*(1-r)*(1-s2))/4.0E+0+(X.xelem(1,9)*(1-r)*(1-s2))/4.0E+0-(X.xelem(1,12)*(1-r2)*(s+1))/4.0E+0+(X.xelem(1,8)*(1-r2)*(s+1))/4.0E+0-(X.xelem(1,14)*(1-r2)*(1-s))/4.0E+0+(X.xelem(1,10)*(1-r2)*(1-s))/4.0E+0;
          J.xelem(0,2) = (-(X.xelem(2,17)*(s+1)*(1-t2))/4.0E+0)+(X.xelem(2,16)*(s+1)*(1-t2))/4.0E+0+(X.xelem(2,19)*(1-s)*(1-t2))/4.0E+0-(X.xelem(2,18)*(1-s)*(1-t2))/4.0E+0+X.xelem(2,0)*(((s+1)*(t+1))/8.0E+0-(((s+1)*(1-t2))/4.0E+0+((1-s2)*(t+1))/4.0E+0-(r*(s+1)*(t+1))/2.0E+0)/2.0E+0)+X.xelem(2,4)*(((s+1)*(1-t))/8.0E+0-(((s+1)*(1-t2))/4.0E+0+((1-s2)*(1-t))/4.0E+0-(r*(s+1)*(1-t))/2.0E+0)/2.0E+0)+X.xelem(2,1)*((-((-((s+1)*(1-t2))/4.0E+0)-((1-s2)*(t+1))/4.0E+0-(r*(s+1)*(t+1))/2.0E+0)/2.0E+0)-((s+1)*(t+1))/8.0E+0)+X.xelem(2,5)*((-((-((s+1)*(1-t2))/4.0E+0)-((1-s2)*(1-t))/4.0E+0-(r*(s+1)*(1-t))/2.0E+0)/2.0E+0)-((s+1)*(1-t))/8.0E+0)+X.xelem(2,3)*(((1-s)*(t+1))/8.0E+0-(((1-s)*(1-t2))/4.0E+0+((1-s2)*(t+1))/4.0E+0-(r*(1-s)*(t+1))/2.0E+0)/2.0E+0)+X.xelem(2,7)*(((1-s)*(1-t))/8.0E+0-(((1-s)*(1-t2))/4.0E+0+((1-s2)*(1-t))/4.0E+0-(r*(1-s)*(1-t))/2.0E+0)/2.0E+0)+X.xelem(2,2)*((-((-((1-s)*(1-t2))/4.0E+0)-((1-s2)*(t+1))/4.0E+0-(r*(1-s)*(t+1))/2.0E+0)/2.0E+0)-((1-s)*(t+1))/8.0E+0)+X.xelem(2,6)*((-((-((1-s)*(1-t2))/4.0E+0)-((1-s2)*(1-t))/4.0E+0-(r*(1-s)*(1-t))/2.0E+0)/2.0E+0)-((1-s)*(1-t))/8.0E+0)+(X.xelem(2,11)*(1-s2)*(t+1))/4.0E+0-(X.xelem(2,9)*(1-s2)*(t+1))/4.0E+0-(X.xelem(2,8)*r*(s+1)*(t+1))/2.0E+0-(X.xelem(2,10)*r*(1-s)*(t+1))/2.0E+0+(X.xelem(2,15)*(1-s2)*(1-t))/4.0E+0-(X.xelem(2,13)*(1-s2)*(1-t))/4.0E+0-(X.xelem(2,12)*r*(s+1)*(1-t))/2.0E+0-(X.xelem(2,14)*r*(1-s)*(1-t))/2.0E+0;
          J.xelem(1,2) = (-(X.xelem(2,19)*(r+1)*(1-t2))/4.0E+0)+(X.xelem(2,16)*(r+1)*(1-t2))/4.0E+0-(X.xelem(2,18)*(1-r)*(1-t2))/4.0E+0+(X.xelem(2,17)*(1-r)*(1-t2))/4.0E+0+X.xelem(2,0)*(((r+1)*(t+1))/8.0E+0-(((r+1)*(1-t2))/4.0E+0-((r+1)*s*(t+1))/2.0E+0+((1-r2)*(t+1))/4.0E+0)/2.0E+0)+X.xelem(2,4)*(((r+1)*(1-t))/8.0E+0-(((r+1)*(1-t2))/4.0E+0-((r+1)*s*(1-t))/2.0E+0+((1-r2)*(1-t))/4.0E+0)/2.0E+0)+X.xelem(2,3)*((-((-((r+1)*(1-t2))/4.0E+0)-((r+1)*s*(t+1))/2.0E+0-((1-r2)*(t+1))/4.0E+0)/2.0E+0)-((r+1)*(t+1))/8.0E+0)+X.xelem(2,7)*((-((-((r+1)*(1-t2))/4.0E+0)-((r+1)*s*(1-t))/2.0E+0-((1-r2)*(1-t))/4.0E+0)/2.0E+0)-((r+1)*(1-t))/8.0E+0)+X.xelem(2,1)*(((1-r)*(t+1))/8.0E+0-(((1-r)*(1-t2))/4.0E+0-((1-r)*s*(t+1))/2.0E+0+((1-r2)*(t+1))/4.0E+0)/2.0E+0)+X.xelem(2,5)*(((1-r)*(1-t))/8.0E+0-(((1-r)*(1-t2))/4.0E+0-((1-r)*s*(1-t))/2.0E+0+((1-r2)*(1-t))/4.0E+0)/2.0E+0)+X.xelem(2,2)*((-((-((1-r)*(1-t2))/4.0E+0)-((1-r)*s*(t+1))/2.0E+0-((1-r2)*(t+1))/4.0E+0)/2.0E+0)-((1-r)*(t+1))/8.0E+0)+X.xelem(2,6)*((-((-((1-r)*(1-t2))/4.0E+0)-((1-r)*s*(1-t))/2.0E+0-((1-r2)*(1-t))/4.0E+0)/2.0E+0)-((1-r)*(1-t))/8.0E+0)-(X.xelem(2,11)*(r+1)*s*(t+1))/2.0E+0-(X.xelem(2,9)*(1-r)*s*(t+1))/2.0E+0-(X.xelem(2,10)*(1-r2)*(t+1))/4.0E+0+(X.xelem(2,8)*(1-r2)*(t+1))/4.0E+0-(X.xelem(2,15)*(r+1)*s*(1-t))/2.0E+0-(X.xelem(2,13)*(1-r)*s*(1-t))/2.0E+0-(X.xelem(2,14)*(1-r2)*(1-t))/4.0E+0+(X.xelem(2,12)*(1-r2)*(1-t))/4.0E+0;
          J.xelem(2,2) = X.xelem(2,0)*(((r+1)*(s+1))/8.0E+0-((-((r+1)*(s+1)*t)/2.0E+0)+((r+1)*(1-s2))/4.0E+0+((1-r2)*(s+1))/4.0E+0)/2.0E+0)+X.xelem(2,4)*((-((-((r+1)*(s+1)*t)/2.0E+0)-((r+1)*(1-s2))/4.0E+0-((1-r2)*(s+1))/4.0E+0)/2.0E+0)-((r+1)*(s+1))/8.0E+0)+X.xelem(2,1)*(((1-r)*(s+1))/8.0E+0-((-((1-r)*(s+1)*t)/2.0E+0)+((1-r)*(1-s2))/4.0E+0+((1-r2)*(s+1))/4.0E+0)/2.0E+0)+X.xelem(2,5)*((-((-((1-r)*(s+1)*t)/2.0E+0)-((1-r)*(1-s2))/4.0E+0-((1-r2)*(s+1))/4.0E+0)/2.0E+0)-((1-r)*(s+1))/8.0E+0)+X.xelem(2,3)*(((r+1)*(1-s))/8.0E+0-((-((r+1)*(1-s)*t)/2.0E+0)+((r+1)*(1-s2))/4.0E+0+((1-r2)*(1-s))/4.0E+0)/2.0E+0)+X.xelem(2,7)*((-((-((r+1)*(1-s)*t)/2.0E+0)-((r+1)*(1-s2))/4.0E+0-((1-r2)*(1-s))/4.0E+0)/2.0E+0)-((r+1)*(1-s))/8.0E+0)+X.xelem(2,2)*(((1-r)*(1-s))/8.0E+0-((-((1-r)*(1-s)*t)/2.0E+0)+((1-r)*(1-s2))/4.0E+0+((1-r2)*(1-s))/4.0E+0)/2.0E+0)+X.xelem(2,6)*((-((-((1-r)*(1-s)*t)/2.0E+0)-((1-r)*(1-s2))/4.0E+0-((1-r2)*(1-s))/4.0E+0)/2.0E+0)-((1-r)*(1-s))/8.0E+0)-(X.xelem(2,16)*(r+1)*(s+1)*t)/2.0E+0-(X.xelem(2,17)*(1-r)*(s+1)*t)/2.0E+0-(X.xelem(2,19)*(r+1)*(1-s)*t)/2.0E+0-(X.xelem(2,18)*(1-r)*(1-s)*t)/2.0E+0-(X.xelem(2,15)*(r+1)*(1-s2))/4.0E+0+(X.xelem(2,11)*(r+1)*(1-s2))/4.0E+0-(X.xelem(2,13)*(1-r)*(1-s2))/4.0E+0+(X.xelem(2,9)*(1-r)*(1-s2))/4.0E+0-(X.xelem(2,12)*(1-r2)*(s+1))/4.0E+0+(X.xelem(2,8)*(1-r2)*(s+1))/4.0E+0-(X.xelem(2,14)*(1-r2)*(1-s))/4.0E+0+(X.xelem(2,10)*(1-r2)*(1-s))/4.0E+0;

          return J.determinant();
     }

     virtual void DispInterpMatrix(const ColumnVector& rv, Matrix& H) const {
          FEM_ASSERT(rv.numel() == 3);
          FEM_ASSERT(H.rows() == 3);
          FEM_ASSERT(H.columns() == 60);

          const double r = rv.xelem(0);
          const double s = rv.xelem(1);
          const double t = rv.xelem(2);
          const double r2 = r * r;
          const double s2 = s * s;
          const double t2 = t * t;

          H.xelem(0,0) = ((r+1)*(s+1)*(t+1))/8.0E+0-(((r+1)*(s+1)*(1-t2))/4.0E+0+((r+1)*(1-s2)*(t+1))/4.0E+0+((1-r2)*(s+1)*(t+1))/4.0E+0)/2.0E+0;
          H.xelem(1,0) = 0;
          H.xelem(2,0) = 0;
          H.xelem(0,1) = 0;
          H.xelem(1,1) = ((r+1)*(s+1)*(t+1))/8.0E+0-(((r+1)*(s+1)*(1-t2))/4.0E+0+((r+1)*(1-s2)*(t+1))/4.0E+0+((1-r2)*(s+1)*(t+1))/4.0E+0)/2.0E+0;
          H.xelem(2,1) = 0;
          H.xelem(0,2) = 0;
          H.xelem(1,2) = 0;
          H.xelem(2,2) = ((r+1)*(s+1)*(t+1))/8.0E+0-(((r+1)*(s+1)*(1-t2))/4.0E+0+((r+1)*(1-s2)*(t+1))/4.0E+0+((1-r2)*(s+1)*(t+1))/4.0E+0)/2.0E+0;
          H.xelem(0,3) = ((1-r)*(s+1)*(t+1))/8.0E+0-(((1-r)*(s+1)*(1-t2))/4.0E+0+((1-r)*(1-s2)*(t+1))/4.0E+0+((1-r2)*(s+1)*(t+1))/4.0E+0)/2.0E+0;
          H.xelem(1,3) = 0;
          H.xelem(2,3) = 0;
          H.xelem(0,4) = 0;
          H.xelem(1,4) = ((1-r)*(s+1)*(t+1))/8.0E+0-(((1-r)*(s+1)*(1-t2))/4.0E+0+((1-r)*(1-s2)*(t+1))/4.0E+0+((1-r2)*(s+1)*(t+1))/4.0E+0)/2.0E+0;
          H.xelem(2,4) = 0;
          H.xelem(0,5) = 0;
          H.xelem(1,5) = 0;
          H.xelem(2,5) = ((1-r)*(s+1)*(t+1))/8.0E+0-(((1-r)*(s+1)*(1-t2))/4.0E+0+((1-r)*(1-s2)*(t+1))/4.0E+0+((1-r2)*(s+1)*(t+1))/4.0E+0)/2.0E+0;
          H.xelem(0,6) = ((1-r)*(1-s)*(t+1))/8.0E+0-(((1-r)*(1-s)*(1-t2))/4.0E+0+((1-r)*(1-s2)*(t+1))/4.0E+0+((1-r2)*(1-s)*(t+1))/4.0E+0)/2.0E+0;
          H.xelem(1,6) = 0;
          H.xelem(2,6) = 0;
          H.xelem(0,7) = 0;
          H.xelem(1,7) = ((1-r)*(1-s)*(t+1))/8.0E+0-(((1-r)*(1-s)*(1-t2))/4.0E+0+((1-r)*(1-s2)*(t+1))/4.0E+0+((1-r2)*(1-s)*(t+1))/4.0E+0)/2.0E+0;
          H.xelem(2,7) = 0;
          H.xelem(0,8) = 0;
          H.xelem(1,8) = 0;
          H.xelem(2,8) = ((1-r)*(1-s)*(t+1))/8.0E+0-(((1-r)*(1-s)*(1-t2))/4.0E+0+((1-r)*(1-s2)*(t+1))/4.0E+0+((1-r2)*(1-s)*(t+1))/4.0E+0)/2.0E+0;
          H.xelem(0,9) = ((r+1)*(1-s)*(t+1))/8.0E+0-(((r+1)*(1-s)*(1-t2))/4.0E+0+((r+1)*(1-s2)*(t+1))/4.0E+0+((1-r2)*(1-s)*(t+1))/4.0E+0)/2.0E+0;
          H.xelem(1,9) = 0;
          H.xelem(2,9) = 0;
          H.xelem(0,10) = 0;
          H.xelem(1,10) = ((r+1)*(1-s)*(t+1))/8.0E+0-(((r+1)*(1-s)*(1-t2))/4.0E+0+((r+1)*(1-s2)*(t+1))/4.0E+0+((1-r2)*(1-s)*(t+1))/4.0E+0)/2.0E+0;
          H.xelem(2,10) = 0;
          H.xelem(0,11) = 0;
          H.xelem(1,11) = 0;
          H.xelem(2,11) = ((r+1)*(1-s)*(t+1))/8.0E+0-(((r+1)*(1-s)*(1-t2))/4.0E+0+((r+1)*(1-s2)*(t+1))/4.0E+0+((1-r2)*(1-s)*(t+1))/4.0E+0)/2.0E+0;
          H.xelem(0,12) = ((r+1)*(s+1)*(1-t))/8.0E+0-(((r+1)*(s+1)*(1-t2))/4.0E+0+((r+1)*(1-s2)*(1-t))/4.0E+0+((1-r2)*(s+1)*(1-t))/4.0E+0)/2.0E+0;
          H.xelem(1,12) = 0;
          H.xelem(2,12) = 0;
          H.xelem(0,13) = 0;
          H.xelem(1,13) = ((r+1)*(s+1)*(1-t))/8.0E+0-(((r+1)*(s+1)*(1-t2))/4.0E+0+((r+1)*(1-s2)*(1-t))/4.0E+0+((1-r2)*(s+1)*(1-t))/4.0E+0)/2.0E+0;
          H.xelem(2,13) = 0;
          H.xelem(0,14) = 0;
          H.xelem(1,14) = 0;
          H.xelem(2,14) = ((r+1)*(s+1)*(1-t))/8.0E+0-(((r+1)*(s+1)*(1-t2))/4.0E+0+((r+1)*(1-s2)*(1-t))/4.0E+0+((1-r2)*(s+1)*(1-t))/4.0E+0)/2.0E+0;
          H.xelem(0,15) = ((1-r)*(s+1)*(1-t))/8.0E+0-(((1-r)*(s+1)*(1-t2))/4.0E+0+((1-r)*(1-s2)*(1-t))/4.0E+0+((1-r2)*(s+1)*(1-t))/4.0E+0)/2.0E+0;
          H.xelem(1,15) = 0;
          H.xelem(2,15) = 0;
          H.xelem(0,16) = 0;
          H.xelem(1,16) = ((1-r)*(s+1)*(1-t))/8.0E+0-(((1-r)*(s+1)*(1-t2))/4.0E+0+((1-r)*(1-s2)*(1-t))/4.0E+0+((1-r2)*(s+1)*(1-t))/4.0E+0)/2.0E+0;
          H.xelem(2,16) = 0;
          H.xelem(0,17) = 0;
          H.xelem(1,17) = 0;
          H.xelem(2,17) = ((1-r)*(s+1)*(1-t))/8.0E+0-(((1-r)*(s+1)*(1-t2))/4.0E+0+((1-r)*(1-s2)*(1-t))/4.0E+0+((1-r2)*(s+1)*(1-t))/4.0E+0)/2.0E+0;
          H.xelem(0,18) = ((1-r)*(1-s)*(1-t))/8.0E+0-(((1-r)*(1-s)*(1-t2))/4.0E+0+((1-r)*(1-s2)*(1-t))/4.0E+0+((1-r2)*(1-s)*(1-t))/4.0E+0)/2.0E+0;
          H.xelem(1,18) = 0;
          H.xelem(2,18) = 0;
          H.xelem(0,19) = 0;
          H.xelem(1,19) = ((1-r)*(1-s)*(1-t))/8.0E+0-(((1-r)*(1-s)*(1-t2))/4.0E+0+((1-r)*(1-s2)*(1-t))/4.0E+0+((1-r2)*(1-s)*(1-t))/4.0E+0)/2.0E+0;
          H.xelem(2,19) = 0;
          H.xelem(0,20) = 0;
          H.xelem(1,20) = 0;
          H.xelem(2,20) = ((1-r)*(1-s)*(1-t))/8.0E+0-(((1-r)*(1-s)*(1-t2))/4.0E+0+((1-r)*(1-s2)*(1-t))/4.0E+0+((1-r2)*(1-s)*(1-t))/4.0E+0)/2.0E+0;
          H.xelem(0,21) = ((r+1)*(1-s)*(1-t))/8.0E+0-(((r+1)*(1-s)*(1-t2))/4.0E+0+((r+1)*(1-s2)*(1-t))/4.0E+0+((1-r2)*(1-s)*(1-t))/4.0E+0)/2.0E+0;
          H.xelem(1,21) = 0;
          H.xelem(2,21) = 0;
          H.xelem(0,22) = 0;
          H.xelem(1,22) = ((r+1)*(1-s)*(1-t))/8.0E+0-(((r+1)*(1-s)*(1-t2))/4.0E+0+((r+1)*(1-s2)*(1-t))/4.0E+0+((1-r2)*(1-s)*(1-t))/4.0E+0)/2.0E+0;
          H.xelem(2,22) = 0;
          H.xelem(0,23) = 0;
          H.xelem(1,23) = 0;
          H.xelem(2,23) = ((r+1)*(1-s)*(1-t))/8.0E+0-(((r+1)*(1-s)*(1-t2))/4.0E+0+((r+1)*(1-s2)*(1-t))/4.0E+0+((1-r2)*(1-s)*(1-t))/4.0E+0)/2.0E+0;
          H.xelem(0,24) = ((1-r2)*(s+1)*(t+1))/4.0E+0;
          H.xelem(1,24) = 0;
          H.xelem(2,24) = 0;
          H.xelem(0,25) = 0;
          H.xelem(1,25) = ((1-r2)*(s+1)*(t+1))/4.0E+0;
          H.xelem(2,25) = 0;
          H.xelem(0,26) = 0;
          H.xelem(1,26) = 0;
          H.xelem(2,26) = ((1-r2)*(s+1)*(t+1))/4.0E+0;
          H.xelem(0,27) = ((1-r)*(1-s2)*(t+1))/4.0E+0;
          H.xelem(1,27) = 0;
          H.xelem(2,27) = 0;
          H.xelem(0,28) = 0;
          H.xelem(1,28) = ((1-r)*(1-s2)*(t+1))/4.0E+0;
          H.xelem(2,28) = 0;
          H.xelem(0,29) = 0;
          H.xelem(1,29) = 0;
          H.xelem(2,29) = ((1-r)*(1-s2)*(t+1))/4.0E+0;
          H.xelem(0,30) = ((1-r2)*(1-s)*(t+1))/4.0E+0;
          H.xelem(1,30) = 0;
          H.xelem(2,30) = 0;
          H.xelem(0,31) = 0;
          H.xelem(1,31) = ((1-r2)*(1-s)*(t+1))/4.0E+0;
          H.xelem(2,31) = 0;
          H.xelem(0,32) = 0;
          H.xelem(1,32) = 0;
          H.xelem(2,32) = ((1-r2)*(1-s)*(t+1))/4.0E+0;
          H.xelem(0,33) = ((r+1)*(1-s2)*(t+1))/4.0E+0;
          H.xelem(1,33) = 0;
          H.xelem(2,33) = 0;
          H.xelem(0,34) = 0;
          H.xelem(1,34) = ((r+1)*(1-s2)*(t+1))/4.0E+0;
          H.xelem(2,34) = 0;
          H.xelem(0,35) = 0;
          H.xelem(1,35) = 0;
          H.xelem(2,35) = ((r+1)*(1-s2)*(t+1))/4.0E+0;
          H.xelem(0,36) = ((1-r2)*(s+1)*(1-t))/4.0E+0;
          H.xelem(1,36) = 0;
          H.xelem(2,36) = 0;
          H.xelem(0,37) = 0;
          H.xelem(1,37) = ((1-r2)*(s+1)*(1-t))/4.0E+0;
          H.xelem(2,37) = 0;
          H.xelem(0,38) = 0;
          H.xelem(1,38) = 0;
          H.xelem(2,38) = ((1-r2)*(s+1)*(1-t))/4.0E+0;
          H.xelem(0,39) = ((1-r)*(1-s2)*(1-t))/4.0E+0;
          H.xelem(1,39) = 0;
          H.xelem(2,39) = 0;
          H.xelem(0,40) = 0;
          H.xelem(1,40) = ((1-r)*(1-s2)*(1-t))/4.0E+0;
          H.xelem(2,40) = 0;
          H.xelem(0,41) = 0;
          H.xelem(1,41) = 0;
          H.xelem(2,41) = ((1-r)*(1-s2)*(1-t))/4.0E+0;
          H.xelem(0,42) = ((1-r2)*(1-s)*(1-t))/4.0E+0;
          H.xelem(1,42) = 0;
          H.xelem(2,42) = 0;
          H.xelem(0,43) = 0;
          H.xelem(1,43) = ((1-r2)*(1-s)*(1-t))/4.0E+0;
          H.xelem(2,43) = 0;
          H.xelem(0,44) = 0;
          H.xelem(1,44) = 0;
          H.xelem(2,44) = ((1-r2)*(1-s)*(1-t))/4.0E+0;
          H.xelem(0,45) = ((r+1)*(1-s2)*(1-t))/4.0E+0;
          H.xelem(1,45) = 0;
          H.xelem(2,45) = 0;
          H.xelem(0,46) = 0;
          H.xelem(1,46) = ((r+1)*(1-s2)*(1-t))/4.0E+0;
          H.xelem(2,46) = 0;
          H.xelem(0,47) = 0;
          H.xelem(1,47) = 0;
          H.xelem(2,47) = ((r+1)*(1-s2)*(1-t))/4.0E+0;
          H.xelem(0,48) = ((r+1)*(s+1)*(1-t2))/4.0E+0;
          H.xelem(1,48) = 0;
          H.xelem(2,48) = 0;
          H.xelem(0,49) = 0;
          H.xelem(1,49) = ((r+1)*(s+1)*(1-t2))/4.0E+0;
          H.xelem(2,49) = 0;
          H.xelem(0,50) = 0;
          H.xelem(1,50) = 0;
          H.xelem(2,50) = ((r+1)*(s+1)*(1-t2))/4.0E+0;
          H.xelem(0,51) = ((1-r)*(s+1)*(1-t2))/4.0E+0;
          H.xelem(1,51) = 0;
          H.xelem(2,51) = 0;
          H.xelem(0,52) = 0;
          H.xelem(1,52) = ((1-r)*(s+1)*(1-t2))/4.0E+0;
          H.xelem(2,52) = 0;
          H.xelem(0,53) = 0;
          H.xelem(1,53) = 0;
          H.xelem(2,53) = ((1-r)*(s+1)*(1-t2))/4.0E+0;
          H.xelem(0,54) = ((1-r)*(1-s)*(1-t2))/4.0E+0;
          H.xelem(1,54) = 0;
          H.xelem(2,54) = 0;
          H.xelem(0,55) = 0;
          H.xelem(1,55) = ((1-r)*(1-s)*(1-t2))/4.0E+0;
          H.xelem(2,55) = 0;
          H.xelem(0,56) = 0;
          H.xelem(1,56) = 0;
          H.xelem(2,56) = ((1-r)*(1-s)*(1-t2))/4.0E+0;
          H.xelem(0,57) = ((r+1)*(1-s)*(1-t2))/4.0E+0;
          H.xelem(1,57) = 0;
          H.xelem(2,57) = 0;
          H.xelem(0,58) = 0;
          H.xelem(1,58) = ((r+1)*(1-s)*(1-t2))/4.0E+0;
          H.xelem(2,58) = 0;
          H.xelem(0,59) = 0;
          H.xelem(1,59) = 0;
          H.xelem(2,59) = ((r+1)*(1-s)*(1-t2))/4.0E+0;
     }

     virtual void StrainMatrix(const ColumnVector& rv, const Matrix& invJ, Matrix& B) const {
          FEM_ASSERT(rv.numel() == 3);
          FEM_ASSERT(invJ.rows() == 3);
          FEM_ASSERT(invJ.columns() == 3);
          FEM_ASSERT(B.rows() == 6);
          FEM_ASSERT(B.columns() == 60);

          const double r = rv.xelem(0);
          const double s = rv.xelem(1);
          const double t = rv.xelem(2);
          const double r2 = r * r;
          const double s2 = s * s;
          const double t2 = t * t;

          B.xelem(0,0) = invJ.xelem(0,0)*(((s+1)*(t+1))/8.0E+0-(((s+1)*(1-t2))/4.0E+0+((1-s2)*(t+1))/4.0E+0-(r*(s+1)*(t+1))/2.0E+0)/2.0E+0)+invJ.xelem(0,1)*(((r+1)*(t+1))/8.0E+0-(((r+1)*(1-t2))/4.0E+0-((r+1)*s*(t+1))/2.0E+0+((1-r2)*(t+1))/4.0E+0)/2.0E+0)+invJ.xelem(0,2)*(((r+1)*(s+1))/8.0E+0-((-((r+1)*(s+1)*t)/2.0E+0)+((r+1)*(1-s2))/4.0E+0+((1-r2)*(s+1))/4.0E+0)/2.0E+0);
          B.xelem(1,0) = 0;
          B.xelem(2,0) = 0;
          B.xelem(3,0) = invJ.xelem(1,0)*(((s+1)*(t+1))/8.0E+0-(((s+1)*(1-t2))/4.0E+0+((1-s2)*(t+1))/4.0E+0-(r*(s+1)*(t+1))/2.0E+0)/2.0E+0)+invJ.xelem(1,1)*(((r+1)*(t+1))/8.0E+0-(((r+1)*(1-t2))/4.0E+0-((r+1)*s*(t+1))/2.0E+0+((1-r2)*(t+1))/4.0E+0)/2.0E+0)+invJ.xelem(1,2)*(((r+1)*(s+1))/8.0E+0-((-((r+1)*(s+1)*t)/2.0E+0)+((r+1)*(1-s2))/4.0E+0+((1-r2)*(s+1))/4.0E+0)/2.0E+0);
          B.xelem(4,0) = 0;
          B.xelem(5,0) = invJ.xelem(2,0)*(((s+1)*(t+1))/8.0E+0-(((s+1)*(1-t2))/4.0E+0+((1-s2)*(t+1))/4.0E+0-(r*(s+1)*(t+1))/2.0E+0)/2.0E+0)+invJ.xelem(2,1)*(((r+1)*(t+1))/8.0E+0-(((r+1)*(1-t2))/4.0E+0-((r+1)*s*(t+1))/2.0E+0+((1-r2)*(t+1))/4.0E+0)/2.0E+0)+invJ.xelem(2,2)*(((r+1)*(s+1))/8.0E+0-((-((r+1)*(s+1)*t)/2.0E+0)+((r+1)*(1-s2))/4.0E+0+((1-r2)*(s+1))/4.0E+0)/2.0E+0);
          B.xelem(0,1) = 0;
          B.xelem(1,1) = invJ.xelem(1,0)*(((s+1)*(t+1))/8.0E+0-(((s+1)*(1-t2))/4.0E+0+((1-s2)*(t+1))/4.0E+0-(r*(s+1)*(t+1))/2.0E+0)/2.0E+0)+invJ.xelem(1,1)*(((r+1)*(t+1))/8.0E+0-(((r+1)*(1-t2))/4.0E+0-((r+1)*s*(t+1))/2.0E+0+((1-r2)*(t+1))/4.0E+0)/2.0E+0)+invJ.xelem(1,2)*(((r+1)*(s+1))/8.0E+0-((-((r+1)*(s+1)*t)/2.0E+0)+((r+1)*(1-s2))/4.0E+0+((1-r2)*(s+1))/4.0E+0)/2.0E+0);
          B.xelem(2,1) = 0;
          B.xelem(3,1) = invJ.xelem(0,0)*(((s+1)*(t+1))/8.0E+0-(((s+1)*(1-t2))/4.0E+0+((1-s2)*(t+1))/4.0E+0-(r*(s+1)*(t+1))/2.0E+0)/2.0E+0)+invJ.xelem(0,1)*(((r+1)*(t+1))/8.0E+0-(((r+1)*(1-t2))/4.0E+0-((r+1)*s*(t+1))/2.0E+0+((1-r2)*(t+1))/4.0E+0)/2.0E+0)+invJ.xelem(0,2)*(((r+1)*(s+1))/8.0E+0-((-((r+1)*(s+1)*t)/2.0E+0)+((r+1)*(1-s2))/4.0E+0+((1-r2)*(s+1))/4.0E+0)/2.0E+0);
          B.xelem(4,1) = invJ.xelem(2,0)*(((s+1)*(t+1))/8.0E+0-(((s+1)*(1-t2))/4.0E+0+((1-s2)*(t+1))/4.0E+0-(r*(s+1)*(t+1))/2.0E+0)/2.0E+0)+invJ.xelem(2,1)*(((r+1)*(t+1))/8.0E+0-(((r+1)*(1-t2))/4.0E+0-((r+1)*s*(t+1))/2.0E+0+((1-r2)*(t+1))/4.0E+0)/2.0E+0)+invJ.xelem(2,2)*(((r+1)*(s+1))/8.0E+0-((-((r+1)*(s+1)*t)/2.0E+0)+((r+1)*(1-s2))/4.0E+0+((1-r2)*(s+1))/4.0E+0)/2.0E+0);
          B.xelem(5,1) = 0;
          B.xelem(0,2) = 0;
          B.xelem(1,2) = 0;
          B.xelem(2,2) = invJ.xelem(2,0)*(((s+1)*(t+1))/8.0E+0-(((s+1)*(1-t2))/4.0E+0+((1-s2)*(t+1))/4.0E+0-(r*(s+1)*(t+1))/2.0E+0)/2.0E+0)+invJ.xelem(2,1)*(((r+1)*(t+1))/8.0E+0-(((r+1)*(1-t2))/4.0E+0-((r+1)*s*(t+1))/2.0E+0+((1-r2)*(t+1))/4.0E+0)/2.0E+0)+invJ.xelem(2,2)*(((r+1)*(s+1))/8.0E+0-((-((r+1)*(s+1)*t)/2.0E+0)+((r+1)*(1-s2))/4.0E+0+((1-r2)*(s+1))/4.0E+0)/2.0E+0);
          B.xelem(3,2) = 0;
          B.xelem(4,2) = invJ.xelem(1,0)*(((s+1)*(t+1))/8.0E+0-(((s+1)*(1-t2))/4.0E+0+((1-s2)*(t+1))/4.0E+0-(r*(s+1)*(t+1))/2.0E+0)/2.0E+0)+invJ.xelem(1,1)*(((r+1)*(t+1))/8.0E+0-(((r+1)*(1-t2))/4.0E+0-((r+1)*s*(t+1))/2.0E+0+((1-r2)*(t+1))/4.0E+0)/2.0E+0)+invJ.xelem(1,2)*(((r+1)*(s+1))/8.0E+0-((-((r+1)*(s+1)*t)/2.0E+0)+((r+1)*(1-s2))/4.0E+0+((1-r2)*(s+1))/4.0E+0)/2.0E+0);
          B.xelem(5,2) = invJ.xelem(0,0)*(((s+1)*(t+1))/8.0E+0-(((s+1)*(1-t2))/4.0E+0+((1-s2)*(t+1))/4.0E+0-(r*(s+1)*(t+1))/2.0E+0)/2.0E+0)+invJ.xelem(0,1)*(((r+1)*(t+1))/8.0E+0-(((r+1)*(1-t2))/4.0E+0-((r+1)*s*(t+1))/2.0E+0+((1-r2)*(t+1))/4.0E+0)/2.0E+0)+invJ.xelem(0,2)*(((r+1)*(s+1))/8.0E+0-((-((r+1)*(s+1)*t)/2.0E+0)+((r+1)*(1-s2))/4.0E+0+((1-r2)*(s+1))/4.0E+0)/2.0E+0);
          B.xelem(0,3) = invJ.xelem(0,0)*((-((-((s+1)*(1-t2))/4.0E+0)-((1-s2)*(t+1))/4.0E+0-(r*(s+1)*(t+1))/2.0E+0)/2.0E+0)-((s+1)*(t+1))/8.0E+0)+invJ.xelem(0,1)*(((1-r)*(t+1))/8.0E+0-(((1-r)*(1-t2))/4.0E+0-((1-r)*s*(t+1))/2.0E+0+((1-r2)*(t+1))/4.0E+0)/2.0E+0)+invJ.xelem(0,2)*(((1-r)*(s+1))/8.0E+0-((-((1-r)*(s+1)*t)/2.0E+0)+((1-r)*(1-s2))/4.0E+0+((1-r2)*(s+1))/4.0E+0)/2.0E+0);
          B.xelem(1,3) = 0;
          B.xelem(2,3) = 0;
          B.xelem(3,3) = invJ.xelem(1,0)*((-((-((s+1)*(1-t2))/4.0E+0)-((1-s2)*(t+1))/4.0E+0-(r*(s+1)*(t+1))/2.0E+0)/2.0E+0)-((s+1)*(t+1))/8.0E+0)+invJ.xelem(1,1)*(((1-r)*(t+1))/8.0E+0-(((1-r)*(1-t2))/4.0E+0-((1-r)*s*(t+1))/2.0E+0+((1-r2)*(t+1))/4.0E+0)/2.0E+0)+invJ.xelem(1,2)*(((1-r)*(s+1))/8.0E+0-((-((1-r)*(s+1)*t)/2.0E+0)+((1-r)*(1-s2))/4.0E+0+((1-r2)*(s+1))/4.0E+0)/2.0E+0);
          B.xelem(4,3) = 0;
          B.xelem(5,3) = invJ.xelem(2,0)*((-((-((s+1)*(1-t2))/4.0E+0)-((1-s2)*(t+1))/4.0E+0-(r*(s+1)*(t+1))/2.0E+0)/2.0E+0)-((s+1)*(t+1))/8.0E+0)+invJ.xelem(2,1)*(((1-r)*(t+1))/8.0E+0-(((1-r)*(1-t2))/4.0E+0-((1-r)*s*(t+1))/2.0E+0+((1-r2)*(t+1))/4.0E+0)/2.0E+0)+invJ.xelem(2,2)*(((1-r)*(s+1))/8.0E+0-((-((1-r)*(s+1)*t)/2.0E+0)+((1-r)*(1-s2))/4.0E+0+((1-r2)*(s+1))/4.0E+0)/2.0E+0);
          B.xelem(0,4) = 0;
          B.xelem(1,4) = invJ.xelem(1,0)*((-((-((s+1)*(1-t2))/4.0E+0)-((1-s2)*(t+1))/4.0E+0-(r*(s+1)*(t+1))/2.0E+0)/2.0E+0)-((s+1)*(t+1))/8.0E+0)+invJ.xelem(1,1)*(((1-r)*(t+1))/8.0E+0-(((1-r)*(1-t2))/4.0E+0-((1-r)*s*(t+1))/2.0E+0+((1-r2)*(t+1))/4.0E+0)/2.0E+0)+invJ.xelem(1,2)*(((1-r)*(s+1))/8.0E+0-((-((1-r)*(s+1)*t)/2.0E+0)+((1-r)*(1-s2))/4.0E+0+((1-r2)*(s+1))/4.0E+0)/2.0E+0);
          B.xelem(2,4) = 0;
          B.xelem(3,4) = invJ.xelem(0,0)*((-((-((s+1)*(1-t2))/4.0E+0)-((1-s2)*(t+1))/4.0E+0-(r*(s+1)*(t+1))/2.0E+0)/2.0E+0)-((s+1)*(t+1))/8.0E+0)+invJ.xelem(0,1)*(((1-r)*(t+1))/8.0E+0-(((1-r)*(1-t2))/4.0E+0-((1-r)*s*(t+1))/2.0E+0+((1-r2)*(t+1))/4.0E+0)/2.0E+0)+invJ.xelem(0,2)*(((1-r)*(s+1))/8.0E+0-((-((1-r)*(s+1)*t)/2.0E+0)+((1-r)*(1-s2))/4.0E+0+((1-r2)*(s+1))/4.0E+0)/2.0E+0);
          B.xelem(4,4) = invJ.xelem(2,0)*((-((-((s+1)*(1-t2))/4.0E+0)-((1-s2)*(t+1))/4.0E+0-(r*(s+1)*(t+1))/2.0E+0)/2.0E+0)-((s+1)*(t+1))/8.0E+0)+invJ.xelem(2,1)*(((1-r)*(t+1))/8.0E+0-(((1-r)*(1-t2))/4.0E+0-((1-r)*s*(t+1))/2.0E+0+((1-r2)*(t+1))/4.0E+0)/2.0E+0)+invJ.xelem(2,2)*(((1-r)*(s+1))/8.0E+0-((-((1-r)*(s+1)*t)/2.0E+0)+((1-r)*(1-s2))/4.0E+0+((1-r2)*(s+1))/4.0E+0)/2.0E+0);
          B.xelem(5,4) = 0;
          B.xelem(0,5) = 0;
          B.xelem(1,5) = 0;
          B.xelem(2,5) = invJ.xelem(2,0)*((-((-((s+1)*(1-t2))/4.0E+0)-((1-s2)*(t+1))/4.0E+0-(r*(s+1)*(t+1))/2.0E+0)/2.0E+0)-((s+1)*(t+1))/8.0E+0)+invJ.xelem(2,1)*(((1-r)*(t+1))/8.0E+0-(((1-r)*(1-t2))/4.0E+0-((1-r)*s*(t+1))/2.0E+0+((1-r2)*(t+1))/4.0E+0)/2.0E+0)+invJ.xelem(2,2)*(((1-r)*(s+1))/8.0E+0-((-((1-r)*(s+1)*t)/2.0E+0)+((1-r)*(1-s2))/4.0E+0+((1-r2)*(s+1))/4.0E+0)/2.0E+0);
          B.xelem(3,5) = 0;
          B.xelem(4,5) = invJ.xelem(1,0)*((-((-((s+1)*(1-t2))/4.0E+0)-((1-s2)*(t+1))/4.0E+0-(r*(s+1)*(t+1))/2.0E+0)/2.0E+0)-((s+1)*(t+1))/8.0E+0)+invJ.xelem(1,1)*(((1-r)*(t+1))/8.0E+0-(((1-r)*(1-t2))/4.0E+0-((1-r)*s*(t+1))/2.0E+0+((1-r2)*(t+1))/4.0E+0)/2.0E+0)+invJ.xelem(1,2)*(((1-r)*(s+1))/8.0E+0-((-((1-r)*(s+1)*t)/2.0E+0)+((1-r)*(1-s2))/4.0E+0+((1-r2)*(s+1))/4.0E+0)/2.0E+0);
          B.xelem(5,5) = invJ.xelem(0,0)*((-((-((s+1)*(1-t2))/4.0E+0)-((1-s2)*(t+1))/4.0E+0-(r*(s+1)*(t+1))/2.0E+0)/2.0E+0)-((s+1)*(t+1))/8.0E+0)+invJ.xelem(0,1)*(((1-r)*(t+1))/8.0E+0-(((1-r)*(1-t2))/4.0E+0-((1-r)*s*(t+1))/2.0E+0+((1-r2)*(t+1))/4.0E+0)/2.0E+0)+invJ.xelem(0,2)*(((1-r)*(s+1))/8.0E+0-((-((1-r)*(s+1)*t)/2.0E+0)+((1-r)*(1-s2))/4.0E+0+((1-r2)*(s+1))/4.0E+0)/2.0E+0);
          B.xelem(0,6) = invJ.xelem(0,0)*((-((-((1-s)*(1-t2))/4.0E+0)-((1-s2)*(t+1))/4.0E+0-(r*(1-s)*(t+1))/2.0E+0)/2.0E+0)-((1-s)*(t+1))/8.0E+0)+invJ.xelem(0,1)*((-((-((1-r)*(1-t2))/4.0E+0)-((1-r)*s*(t+1))/2.0E+0-((1-r2)*(t+1))/4.0E+0)/2.0E+0)-((1-r)*(t+1))/8.0E+0)+invJ.xelem(0,2)*(((1-r)*(1-s))/8.0E+0-((-((1-r)*(1-s)*t)/2.0E+0)+((1-r)*(1-s2))/4.0E+0+((1-r2)*(1-s))/4.0E+0)/2.0E+0);
          B.xelem(1,6) = 0;
          B.xelem(2,6) = 0;
          B.xelem(3,6) = invJ.xelem(1,0)*((-((-((1-s)*(1-t2))/4.0E+0)-((1-s2)*(t+1))/4.0E+0-(r*(1-s)*(t+1))/2.0E+0)/2.0E+0)-((1-s)*(t+1))/8.0E+0)+invJ.xelem(1,1)*((-((-((1-r)*(1-t2))/4.0E+0)-((1-r)*s*(t+1))/2.0E+0-((1-r2)*(t+1))/4.0E+0)/2.0E+0)-((1-r)*(t+1))/8.0E+0)+invJ.xelem(1,2)*(((1-r)*(1-s))/8.0E+0-((-((1-r)*(1-s)*t)/2.0E+0)+((1-r)*(1-s2))/4.0E+0+((1-r2)*(1-s))/4.0E+0)/2.0E+0);
          B.xelem(4,6) = 0;
          B.xelem(5,6) = invJ.xelem(2,0)*((-((-((1-s)*(1-t2))/4.0E+0)-((1-s2)*(t+1))/4.0E+0-(r*(1-s)*(t+1))/2.0E+0)/2.0E+0)-((1-s)*(t+1))/8.0E+0)+invJ.xelem(2,1)*((-((-((1-r)*(1-t2))/4.0E+0)-((1-r)*s*(t+1))/2.0E+0-((1-r2)*(t+1))/4.0E+0)/2.0E+0)-((1-r)*(t+1))/8.0E+0)+invJ.xelem(2,2)*(((1-r)*(1-s))/8.0E+0-((-((1-r)*(1-s)*t)/2.0E+0)+((1-r)*(1-s2))/4.0E+0+((1-r2)*(1-s))/4.0E+0)/2.0E+0);
          B.xelem(0,7) = 0;
          B.xelem(1,7) = invJ.xelem(1,0)*((-((-((1-s)*(1-t2))/4.0E+0)-((1-s2)*(t+1))/4.0E+0-(r*(1-s)*(t+1))/2.0E+0)/2.0E+0)-((1-s)*(t+1))/8.0E+0)+invJ.xelem(1,1)*((-((-((1-r)*(1-t2))/4.0E+0)-((1-r)*s*(t+1))/2.0E+0-((1-r2)*(t+1))/4.0E+0)/2.0E+0)-((1-r)*(t+1))/8.0E+0)+invJ.xelem(1,2)*(((1-r)*(1-s))/8.0E+0-((-((1-r)*(1-s)*t)/2.0E+0)+((1-r)*(1-s2))/4.0E+0+((1-r2)*(1-s))/4.0E+0)/2.0E+0);
          B.xelem(2,7) = 0;
          B.xelem(3,7) = invJ.xelem(0,0)*((-((-((1-s)*(1-t2))/4.0E+0)-((1-s2)*(t+1))/4.0E+0-(r*(1-s)*(t+1))/2.0E+0)/2.0E+0)-((1-s)*(t+1))/8.0E+0)+invJ.xelem(0,1)*((-((-((1-r)*(1-t2))/4.0E+0)-((1-r)*s*(t+1))/2.0E+0-((1-r2)*(t+1))/4.0E+0)/2.0E+0)-((1-r)*(t+1))/8.0E+0)+invJ.xelem(0,2)*(((1-r)*(1-s))/8.0E+0-((-((1-r)*(1-s)*t)/2.0E+0)+((1-r)*(1-s2))/4.0E+0+((1-r2)*(1-s))/4.0E+0)/2.0E+0);
          B.xelem(4,7) = invJ.xelem(2,0)*((-((-((1-s)*(1-t2))/4.0E+0)-((1-s2)*(t+1))/4.0E+0-(r*(1-s)*(t+1))/2.0E+0)/2.0E+0)-((1-s)*(t+1))/8.0E+0)+invJ.xelem(2,1)*((-((-((1-r)*(1-t2))/4.0E+0)-((1-r)*s*(t+1))/2.0E+0-((1-r2)*(t+1))/4.0E+0)/2.0E+0)-((1-r)*(t+1))/8.0E+0)+invJ.xelem(2,2)*(((1-r)*(1-s))/8.0E+0-((-((1-r)*(1-s)*t)/2.0E+0)+((1-r)*(1-s2))/4.0E+0+((1-r2)*(1-s))/4.0E+0)/2.0E+0);
          B.xelem(5,7) = 0;
          B.xelem(0,8) = 0;
          B.xelem(1,8) = 0;
          B.xelem(2,8) = invJ.xelem(2,0)*((-((-((1-s)*(1-t2))/4.0E+0)-((1-s2)*(t+1))/4.0E+0-(r*(1-s)*(t+1))/2.0E+0)/2.0E+0)-((1-s)*(t+1))/8.0E+0)+invJ.xelem(2,1)*((-((-((1-r)*(1-t2))/4.0E+0)-((1-r)*s*(t+1))/2.0E+0-((1-r2)*(t+1))/4.0E+0)/2.0E+0)-((1-r)*(t+1))/8.0E+0)+invJ.xelem(2,2)*(((1-r)*(1-s))/8.0E+0-((-((1-r)*(1-s)*t)/2.0E+0)+((1-r)*(1-s2))/4.0E+0+((1-r2)*(1-s))/4.0E+0)/2.0E+0);
          B.xelem(3,8) = 0;
          B.xelem(4,8) = invJ.xelem(1,0)*((-((-((1-s)*(1-t2))/4.0E+0)-((1-s2)*(t+1))/4.0E+0-(r*(1-s)*(t+1))/2.0E+0)/2.0E+0)-((1-s)*(t+1))/8.0E+0)+invJ.xelem(1,1)*((-((-((1-r)*(1-t2))/4.0E+0)-((1-r)*s*(t+1))/2.0E+0-((1-r2)*(t+1))/4.0E+0)/2.0E+0)-((1-r)*(t+1))/8.0E+0)+invJ.xelem(1,2)*(((1-r)*(1-s))/8.0E+0-((-((1-r)*(1-s)*t)/2.0E+0)+((1-r)*(1-s2))/4.0E+0+((1-r2)*(1-s))/4.0E+0)/2.0E+0);
          B.xelem(5,8) = invJ.xelem(0,0)*((-((-((1-s)*(1-t2))/4.0E+0)-((1-s2)*(t+1))/4.0E+0-(r*(1-s)*(t+1))/2.0E+0)/2.0E+0)-((1-s)*(t+1))/8.0E+0)+invJ.xelem(0,1)*((-((-((1-r)*(1-t2))/4.0E+0)-((1-r)*s*(t+1))/2.0E+0-((1-r2)*(t+1))/4.0E+0)/2.0E+0)-((1-r)*(t+1))/8.0E+0)+invJ.xelem(0,2)*(((1-r)*(1-s))/8.0E+0-((-((1-r)*(1-s)*t)/2.0E+0)+((1-r)*(1-s2))/4.0E+0+((1-r2)*(1-s))/4.0E+0)/2.0E+0);
          B.xelem(0,9) = invJ.xelem(0,0)*(((1-s)*(t+1))/8.0E+0-(((1-s)*(1-t2))/4.0E+0+((1-s2)*(t+1))/4.0E+0-(r*(1-s)*(t+1))/2.0E+0)/2.0E+0)+invJ.xelem(0,1)*((-((-((r+1)*(1-t2))/4.0E+0)-((r+1)*s*(t+1))/2.0E+0-((1-r2)*(t+1))/4.0E+0)/2.0E+0)-((r+1)*(t+1))/8.0E+0)+invJ.xelem(0,2)*(((r+1)*(1-s))/8.0E+0-((-((r+1)*(1-s)*t)/2.0E+0)+((r+1)*(1-s2))/4.0E+0+((1-r2)*(1-s))/4.0E+0)/2.0E+0);
          B.xelem(1,9) = 0;
          B.xelem(2,9) = 0;
          B.xelem(3,9) = invJ.xelem(1,0)*(((1-s)*(t+1))/8.0E+0-(((1-s)*(1-t2))/4.0E+0+((1-s2)*(t+1))/4.0E+0-(r*(1-s)*(t+1))/2.0E+0)/2.0E+0)+invJ.xelem(1,1)*((-((-((r+1)*(1-t2))/4.0E+0)-((r+1)*s*(t+1))/2.0E+0-((1-r2)*(t+1))/4.0E+0)/2.0E+0)-((r+1)*(t+1))/8.0E+0)+invJ.xelem(1,2)*(((r+1)*(1-s))/8.0E+0-((-((r+1)*(1-s)*t)/2.0E+0)+((r+1)*(1-s2))/4.0E+0+((1-r2)*(1-s))/4.0E+0)/2.0E+0);
          B.xelem(4,9) = 0;
          B.xelem(5,9) = invJ.xelem(2,0)*(((1-s)*(t+1))/8.0E+0-(((1-s)*(1-t2))/4.0E+0+((1-s2)*(t+1))/4.0E+0-(r*(1-s)*(t+1))/2.0E+0)/2.0E+0)+invJ.xelem(2,1)*((-((-((r+1)*(1-t2))/4.0E+0)-((r+1)*s*(t+1))/2.0E+0-((1-r2)*(t+1))/4.0E+0)/2.0E+0)-((r+1)*(t+1))/8.0E+0)+invJ.xelem(2,2)*(((r+1)*(1-s))/8.0E+0-((-((r+1)*(1-s)*t)/2.0E+0)+((r+1)*(1-s2))/4.0E+0+((1-r2)*(1-s))/4.0E+0)/2.0E+0);
          B.xelem(0,10) = 0;
          B.xelem(1,10) = invJ.xelem(1,0)*(((1-s)*(t+1))/8.0E+0-(((1-s)*(1-t2))/4.0E+0+((1-s2)*(t+1))/4.0E+0-(r*(1-s)*(t+1))/2.0E+0)/2.0E+0)+invJ.xelem(1,1)*((-((-((r+1)*(1-t2))/4.0E+0)-((r+1)*s*(t+1))/2.0E+0-((1-r2)*(t+1))/4.0E+0)/2.0E+0)-((r+1)*(t+1))/8.0E+0)+invJ.xelem(1,2)*(((r+1)*(1-s))/8.0E+0-((-((r+1)*(1-s)*t)/2.0E+0)+((r+1)*(1-s2))/4.0E+0+((1-r2)*(1-s))/4.0E+0)/2.0E+0);
          B.xelem(2,10) = 0;
          B.xelem(3,10) = invJ.xelem(0,0)*(((1-s)*(t+1))/8.0E+0-(((1-s)*(1-t2))/4.0E+0+((1-s2)*(t+1))/4.0E+0-(r*(1-s)*(t+1))/2.0E+0)/2.0E+0)+invJ.xelem(0,1)*((-((-((r+1)*(1-t2))/4.0E+0)-((r+1)*s*(t+1))/2.0E+0-((1-r2)*(t+1))/4.0E+0)/2.0E+0)-((r+1)*(t+1))/8.0E+0)+invJ.xelem(0,2)*(((r+1)*(1-s))/8.0E+0-((-((r+1)*(1-s)*t)/2.0E+0)+((r+1)*(1-s2))/4.0E+0+((1-r2)*(1-s))/4.0E+0)/2.0E+0);
          B.xelem(4,10) = invJ.xelem(2,0)*(((1-s)*(t+1))/8.0E+0-(((1-s)*(1-t2))/4.0E+0+((1-s2)*(t+1))/4.0E+0-(r*(1-s)*(t+1))/2.0E+0)/2.0E+0)+invJ.xelem(2,1)*((-((-((r+1)*(1-t2))/4.0E+0)-((r+1)*s*(t+1))/2.0E+0-((1-r2)*(t+1))/4.0E+0)/2.0E+0)-((r+1)*(t+1))/8.0E+0)+invJ.xelem(2,2)*(((r+1)*(1-s))/8.0E+0-((-((r+1)*(1-s)*t)/2.0E+0)+((r+1)*(1-s2))/4.0E+0+((1-r2)*(1-s))/4.0E+0)/2.0E+0);
          B.xelem(5,10) = 0;
          B.xelem(0,11) = 0;
          B.xelem(1,11) = 0;
          B.xelem(2,11) = invJ.xelem(2,0)*(((1-s)*(t+1))/8.0E+0-(((1-s)*(1-t2))/4.0E+0+((1-s2)*(t+1))/4.0E+0-(r*(1-s)*(t+1))/2.0E+0)/2.0E+0)+invJ.xelem(2,1)*((-((-((r+1)*(1-t2))/4.0E+0)-((r+1)*s*(t+1))/2.0E+0-((1-r2)*(t+1))/4.0E+0)/2.0E+0)-((r+1)*(t+1))/8.0E+0)+invJ.xelem(2,2)*(((r+1)*(1-s))/8.0E+0-((-((r+1)*(1-s)*t)/2.0E+0)+((r+1)*(1-s2))/4.0E+0+((1-r2)*(1-s))/4.0E+0)/2.0E+0);
          B.xelem(3,11) = 0;
          B.xelem(4,11) = invJ.xelem(1,0)*(((1-s)*(t+1))/8.0E+0-(((1-s)*(1-t2))/4.0E+0+((1-s2)*(t+1))/4.0E+0-(r*(1-s)*(t+1))/2.0E+0)/2.0E+0)+invJ.xelem(1,1)*((-((-((r+1)*(1-t2))/4.0E+0)-((r+1)*s*(t+1))/2.0E+0-((1-r2)*(t+1))/4.0E+0)/2.0E+0)-((r+1)*(t+1))/8.0E+0)+invJ.xelem(1,2)*(((r+1)*(1-s))/8.0E+0-((-((r+1)*(1-s)*t)/2.0E+0)+((r+1)*(1-s2))/4.0E+0+((1-r2)*(1-s))/4.0E+0)/2.0E+0);
          B.xelem(5,11) = invJ.xelem(0,0)*(((1-s)*(t+1))/8.0E+0-(((1-s)*(1-t2))/4.0E+0+((1-s2)*(t+1))/4.0E+0-(r*(1-s)*(t+1))/2.0E+0)/2.0E+0)+invJ.xelem(0,1)*((-((-((r+1)*(1-t2))/4.0E+0)-((r+1)*s*(t+1))/2.0E+0-((1-r2)*(t+1))/4.0E+0)/2.0E+0)-((r+1)*(t+1))/8.0E+0)+invJ.xelem(0,2)*(((r+1)*(1-s))/8.0E+0-((-((r+1)*(1-s)*t)/2.0E+0)+((r+1)*(1-s2))/4.0E+0+((1-r2)*(1-s))/4.0E+0)/2.0E+0);
          B.xelem(0,12) = invJ.xelem(0,0)*(((s+1)*(1-t))/8.0E+0-(((s+1)*(1-t2))/4.0E+0+((1-s2)*(1-t))/4.0E+0-(r*(s+1)*(1-t))/2.0E+0)/2.0E+0)+invJ.xelem(0,1)*(((r+1)*(1-t))/8.0E+0-(((r+1)*(1-t2))/4.0E+0-((r+1)*s*(1-t))/2.0E+0+((1-r2)*(1-t))/4.0E+0)/2.0E+0)+invJ.xelem(0,2)*((-((-((r+1)*(s+1)*t)/2.0E+0)-((r+1)*(1-s2))/4.0E+0-((1-r2)*(s+1))/4.0E+0)/2.0E+0)-((r+1)*(s+1))/8.0E+0);
          B.xelem(1,12) = 0;
          B.xelem(2,12) = 0;
          B.xelem(3,12) = invJ.xelem(1,0)*(((s+1)*(1-t))/8.0E+0-(((s+1)*(1-t2))/4.0E+0+((1-s2)*(1-t))/4.0E+0-(r*(s+1)*(1-t))/2.0E+0)/2.0E+0)+invJ.xelem(1,1)*(((r+1)*(1-t))/8.0E+0-(((r+1)*(1-t2))/4.0E+0-((r+1)*s*(1-t))/2.0E+0+((1-r2)*(1-t))/4.0E+0)/2.0E+0)+invJ.xelem(1,2)*((-((-((r+1)*(s+1)*t)/2.0E+0)-((r+1)*(1-s2))/4.0E+0-((1-r2)*(s+1))/4.0E+0)/2.0E+0)-((r+1)*(s+1))/8.0E+0);
          B.xelem(4,12) = 0;
          B.xelem(5,12) = invJ.xelem(2,0)*(((s+1)*(1-t))/8.0E+0-(((s+1)*(1-t2))/4.0E+0+((1-s2)*(1-t))/4.0E+0-(r*(s+1)*(1-t))/2.0E+0)/2.0E+0)+invJ.xelem(2,1)*(((r+1)*(1-t))/8.0E+0-(((r+1)*(1-t2))/4.0E+0-((r+1)*s*(1-t))/2.0E+0+((1-r2)*(1-t))/4.0E+0)/2.0E+0)+invJ.xelem(2,2)*((-((-((r+1)*(s+1)*t)/2.0E+0)-((r+1)*(1-s2))/4.0E+0-((1-r2)*(s+1))/4.0E+0)/2.0E+0)-((r+1)*(s+1))/8.0E+0);
          B.xelem(0,13) = 0;
          B.xelem(1,13) = invJ.xelem(1,0)*(((s+1)*(1-t))/8.0E+0-(((s+1)*(1-t2))/4.0E+0+((1-s2)*(1-t))/4.0E+0-(r*(s+1)*(1-t))/2.0E+0)/2.0E+0)+invJ.xelem(1,1)*(((r+1)*(1-t))/8.0E+0-(((r+1)*(1-t2))/4.0E+0-((r+1)*s*(1-t))/2.0E+0+((1-r2)*(1-t))/4.0E+0)/2.0E+0)+invJ.xelem(1,2)*((-((-((r+1)*(s+1)*t)/2.0E+0)-((r+1)*(1-s2))/4.0E+0-((1-r2)*(s+1))/4.0E+0)/2.0E+0)-((r+1)*(s+1))/8.0E+0);
          B.xelem(2,13) = 0;
          B.xelem(3,13) = invJ.xelem(0,0)*(((s+1)*(1-t))/8.0E+0-(((s+1)*(1-t2))/4.0E+0+((1-s2)*(1-t))/4.0E+0-(r*(s+1)*(1-t))/2.0E+0)/2.0E+0)+invJ.xelem(0,1)*(((r+1)*(1-t))/8.0E+0-(((r+1)*(1-t2))/4.0E+0-((r+1)*s*(1-t))/2.0E+0+((1-r2)*(1-t))/4.0E+0)/2.0E+0)+invJ.xelem(0,2)*((-((-((r+1)*(s+1)*t)/2.0E+0)-((r+1)*(1-s2))/4.0E+0-((1-r2)*(s+1))/4.0E+0)/2.0E+0)-((r+1)*(s+1))/8.0E+0);
          B.xelem(4,13) = invJ.xelem(2,0)*(((s+1)*(1-t))/8.0E+0-(((s+1)*(1-t2))/4.0E+0+((1-s2)*(1-t))/4.0E+0-(r*(s+1)*(1-t))/2.0E+0)/2.0E+0)+invJ.xelem(2,1)*(((r+1)*(1-t))/8.0E+0-(((r+1)*(1-t2))/4.0E+0-((r+1)*s*(1-t))/2.0E+0+((1-r2)*(1-t))/4.0E+0)/2.0E+0)+invJ.xelem(2,2)*((-((-((r+1)*(s+1)*t)/2.0E+0)-((r+1)*(1-s2))/4.0E+0-((1-r2)*(s+1))/4.0E+0)/2.0E+0)-((r+1)*(s+1))/8.0E+0);
          B.xelem(5,13) = 0;
          B.xelem(0,14) = 0;
          B.xelem(1,14) = 0;
          B.xelem(2,14) = invJ.xelem(2,0)*(((s+1)*(1-t))/8.0E+0-(((s+1)*(1-t2))/4.0E+0+((1-s2)*(1-t))/4.0E+0-(r*(s+1)*(1-t))/2.0E+0)/2.0E+0)+invJ.xelem(2,1)*(((r+1)*(1-t))/8.0E+0-(((r+1)*(1-t2))/4.0E+0-((r+1)*s*(1-t))/2.0E+0+((1-r2)*(1-t))/4.0E+0)/2.0E+0)+invJ.xelem(2,2)*((-((-((r+1)*(s+1)*t)/2.0E+0)-((r+1)*(1-s2))/4.0E+0-((1-r2)*(s+1))/4.0E+0)/2.0E+0)-((r+1)*(s+1))/8.0E+0);
          B.xelem(3,14) = 0;
          B.xelem(4,14) = invJ.xelem(1,0)*(((s+1)*(1-t))/8.0E+0-(((s+1)*(1-t2))/4.0E+0+((1-s2)*(1-t))/4.0E+0-(r*(s+1)*(1-t))/2.0E+0)/2.0E+0)+invJ.xelem(1,1)*(((r+1)*(1-t))/8.0E+0-(((r+1)*(1-t2))/4.0E+0-((r+1)*s*(1-t))/2.0E+0+((1-r2)*(1-t))/4.0E+0)/2.0E+0)+invJ.xelem(1,2)*((-((-((r+1)*(s+1)*t)/2.0E+0)-((r+1)*(1-s2))/4.0E+0-((1-r2)*(s+1))/4.0E+0)/2.0E+0)-((r+1)*(s+1))/8.0E+0);
          B.xelem(5,14) = invJ.xelem(0,0)*(((s+1)*(1-t))/8.0E+0-(((s+1)*(1-t2))/4.0E+0+((1-s2)*(1-t))/4.0E+0-(r*(s+1)*(1-t))/2.0E+0)/2.0E+0)+invJ.xelem(0,1)*(((r+1)*(1-t))/8.0E+0-(((r+1)*(1-t2))/4.0E+0-((r+1)*s*(1-t))/2.0E+0+((1-r2)*(1-t))/4.0E+0)/2.0E+0)+invJ.xelem(0,2)*((-((-((r+1)*(s+1)*t)/2.0E+0)-((r+1)*(1-s2))/4.0E+0-((1-r2)*(s+1))/4.0E+0)/2.0E+0)-((r+1)*(s+1))/8.0E+0);
          B.xelem(0,15) = invJ.xelem(0,0)*((-((-((s+1)*(1-t2))/4.0E+0)-((1-s2)*(1-t))/4.0E+0-(r*(s+1)*(1-t))/2.0E+0)/2.0E+0)-((s+1)*(1-t))/8.0E+0)+invJ.xelem(0,1)*(((1-r)*(1-t))/8.0E+0-(((1-r)*(1-t2))/4.0E+0-((1-r)*s*(1-t))/2.0E+0+((1-r2)*(1-t))/4.0E+0)/2.0E+0)+invJ.xelem(0,2)*((-((-((1-r)*(s+1)*t)/2.0E+0)-((1-r)*(1-s2))/4.0E+0-((1-r2)*(s+1))/4.0E+0)/2.0E+0)-((1-r)*(s+1))/8.0E+0);
          B.xelem(1,15) = 0;
          B.xelem(2,15) = 0;
          B.xelem(3,15) = invJ.xelem(1,0)*((-((-((s+1)*(1-t2))/4.0E+0)-((1-s2)*(1-t))/4.0E+0-(r*(s+1)*(1-t))/2.0E+0)/2.0E+0)-((s+1)*(1-t))/8.0E+0)+invJ.xelem(1,1)*(((1-r)*(1-t))/8.0E+0-(((1-r)*(1-t2))/4.0E+0-((1-r)*s*(1-t))/2.0E+0+((1-r2)*(1-t))/4.0E+0)/2.0E+0)+invJ.xelem(1,2)*((-((-((1-r)*(s+1)*t)/2.0E+0)-((1-r)*(1-s2))/4.0E+0-((1-r2)*(s+1))/4.0E+0)/2.0E+0)-((1-r)*(s+1))/8.0E+0);
          B.xelem(4,15) = 0;
          B.xelem(5,15) = invJ.xelem(2,0)*((-((-((s+1)*(1-t2))/4.0E+0)-((1-s2)*(1-t))/4.0E+0-(r*(s+1)*(1-t))/2.0E+0)/2.0E+0)-((s+1)*(1-t))/8.0E+0)+invJ.xelem(2,1)*(((1-r)*(1-t))/8.0E+0-(((1-r)*(1-t2))/4.0E+0-((1-r)*s*(1-t))/2.0E+0+((1-r2)*(1-t))/4.0E+0)/2.0E+0)+invJ.xelem(2,2)*((-((-((1-r)*(s+1)*t)/2.0E+0)-((1-r)*(1-s2))/4.0E+0-((1-r2)*(s+1))/4.0E+0)/2.0E+0)-((1-r)*(s+1))/8.0E+0);
          B.xelem(0,16) = 0;
          B.xelem(1,16) = invJ.xelem(1,0)*((-((-((s+1)*(1-t2))/4.0E+0)-((1-s2)*(1-t))/4.0E+0-(r*(s+1)*(1-t))/2.0E+0)/2.0E+0)-((s+1)*(1-t))/8.0E+0)+invJ.xelem(1,1)*(((1-r)*(1-t))/8.0E+0-(((1-r)*(1-t2))/4.0E+0-((1-r)*s*(1-t))/2.0E+0+((1-r2)*(1-t))/4.0E+0)/2.0E+0)+invJ.xelem(1,2)*((-((-((1-r)*(s+1)*t)/2.0E+0)-((1-r)*(1-s2))/4.0E+0-((1-r2)*(s+1))/4.0E+0)/2.0E+0)-((1-r)*(s+1))/8.0E+0);
          B.xelem(2,16) = 0;
          B.xelem(3,16) = invJ.xelem(0,0)*((-((-((s+1)*(1-t2))/4.0E+0)-((1-s2)*(1-t))/4.0E+0-(r*(s+1)*(1-t))/2.0E+0)/2.0E+0)-((s+1)*(1-t))/8.0E+0)+invJ.xelem(0,1)*(((1-r)*(1-t))/8.0E+0-(((1-r)*(1-t2))/4.0E+0-((1-r)*s*(1-t))/2.0E+0+((1-r2)*(1-t))/4.0E+0)/2.0E+0)+invJ.xelem(0,2)*((-((-((1-r)*(s+1)*t)/2.0E+0)-((1-r)*(1-s2))/4.0E+0-((1-r2)*(s+1))/4.0E+0)/2.0E+0)-((1-r)*(s+1))/8.0E+0);
          B.xelem(4,16) = invJ.xelem(2,0)*((-((-((s+1)*(1-t2))/4.0E+0)-((1-s2)*(1-t))/4.0E+0-(r*(s+1)*(1-t))/2.0E+0)/2.0E+0)-((s+1)*(1-t))/8.0E+0)+invJ.xelem(2,1)*(((1-r)*(1-t))/8.0E+0-(((1-r)*(1-t2))/4.0E+0-((1-r)*s*(1-t))/2.0E+0+((1-r2)*(1-t))/4.0E+0)/2.0E+0)+invJ.xelem(2,2)*((-((-((1-r)*(s+1)*t)/2.0E+0)-((1-r)*(1-s2))/4.0E+0-((1-r2)*(s+1))/4.0E+0)/2.0E+0)-((1-r)*(s+1))/8.0E+0);
          B.xelem(5,16) = 0;
          B.xelem(0,17) = 0;
          B.xelem(1,17) = 0;
          B.xelem(2,17) = invJ.xelem(2,0)*((-((-((s+1)*(1-t2))/4.0E+0)-((1-s2)*(1-t))/4.0E+0-(r*(s+1)*(1-t))/2.0E+0)/2.0E+0)-((s+1)*(1-t))/8.0E+0)+invJ.xelem(2,1)*(((1-r)*(1-t))/8.0E+0-(((1-r)*(1-t2))/4.0E+0-((1-r)*s*(1-t))/2.0E+0+((1-r2)*(1-t))/4.0E+0)/2.0E+0)+invJ.xelem(2,2)*((-((-((1-r)*(s+1)*t)/2.0E+0)-((1-r)*(1-s2))/4.0E+0-((1-r2)*(s+1))/4.0E+0)/2.0E+0)-((1-r)*(s+1))/8.0E+0);
          B.xelem(3,17) = 0;
          B.xelem(4,17) = invJ.xelem(1,0)*((-((-((s+1)*(1-t2))/4.0E+0)-((1-s2)*(1-t))/4.0E+0-(r*(s+1)*(1-t))/2.0E+0)/2.0E+0)-((s+1)*(1-t))/8.0E+0)+invJ.xelem(1,1)*(((1-r)*(1-t))/8.0E+0-(((1-r)*(1-t2))/4.0E+0-((1-r)*s*(1-t))/2.0E+0+((1-r2)*(1-t))/4.0E+0)/2.0E+0)+invJ.xelem(1,2)*((-((-((1-r)*(s+1)*t)/2.0E+0)-((1-r)*(1-s2))/4.0E+0-((1-r2)*(s+1))/4.0E+0)/2.0E+0)-((1-r)*(s+1))/8.0E+0);
          B.xelem(5,17) = invJ.xelem(0,0)*((-((-((s+1)*(1-t2))/4.0E+0)-((1-s2)*(1-t))/4.0E+0-(r*(s+1)*(1-t))/2.0E+0)/2.0E+0)-((s+1)*(1-t))/8.0E+0)+invJ.xelem(0,1)*(((1-r)*(1-t))/8.0E+0-(((1-r)*(1-t2))/4.0E+0-((1-r)*s*(1-t))/2.0E+0+((1-r2)*(1-t))/4.0E+0)/2.0E+0)+invJ.xelem(0,2)*((-((-((1-r)*(s+1)*t)/2.0E+0)-((1-r)*(1-s2))/4.0E+0-((1-r2)*(s+1))/4.0E+0)/2.0E+0)-((1-r)*(s+1))/8.0E+0);
          B.xelem(0,18) = invJ.xelem(0,0)*((-((-((1-s)*(1-t2))/4.0E+0)-((1-s2)*(1-t))/4.0E+0-(r*(1-s)*(1-t))/2.0E+0)/2.0E+0)-((1-s)*(1-t))/8.0E+0)+invJ.xelem(0,1)*((-((-((1-r)*(1-t2))/4.0E+0)-((1-r)*s*(1-t))/2.0E+0-((1-r2)*(1-t))/4.0E+0)/2.0E+0)-((1-r)*(1-t))/8.0E+0)+invJ.xelem(0,2)*((-((-((1-r)*(1-s)*t)/2.0E+0)-((1-r)*(1-s2))/4.0E+0-((1-r2)*(1-s))/4.0E+0)/2.0E+0)-((1-r)*(1-s))/8.0E+0);
          B.xelem(1,18) = 0;
          B.xelem(2,18) = 0;
          B.xelem(3,18) = invJ.xelem(1,0)*((-((-((1-s)*(1-t2))/4.0E+0)-((1-s2)*(1-t))/4.0E+0-(r*(1-s)*(1-t))/2.0E+0)/2.0E+0)-((1-s)*(1-t))/8.0E+0)+invJ.xelem(1,1)*((-((-((1-r)*(1-t2))/4.0E+0)-((1-r)*s*(1-t))/2.0E+0-((1-r2)*(1-t))/4.0E+0)/2.0E+0)-((1-r)*(1-t))/8.0E+0)+invJ.xelem(1,2)*((-((-((1-r)*(1-s)*t)/2.0E+0)-((1-r)*(1-s2))/4.0E+0-((1-r2)*(1-s))/4.0E+0)/2.0E+0)-((1-r)*(1-s))/8.0E+0);
          B.xelem(4,18) = 0;
          B.xelem(5,18) = invJ.xelem(2,0)*((-((-((1-s)*(1-t2))/4.0E+0)-((1-s2)*(1-t))/4.0E+0-(r*(1-s)*(1-t))/2.0E+0)/2.0E+0)-((1-s)*(1-t))/8.0E+0)+invJ.xelem(2,1)*((-((-((1-r)*(1-t2))/4.0E+0)-((1-r)*s*(1-t))/2.0E+0-((1-r2)*(1-t))/4.0E+0)/2.0E+0)-((1-r)*(1-t))/8.0E+0)+invJ.xelem(2,2)*((-((-((1-r)*(1-s)*t)/2.0E+0)-((1-r)*(1-s2))/4.0E+0-((1-r2)*(1-s))/4.0E+0)/2.0E+0)-((1-r)*(1-s))/8.0E+0);
          B.xelem(0,19) = 0;
          B.xelem(1,19) = invJ.xelem(1,0)*((-((-((1-s)*(1-t2))/4.0E+0)-((1-s2)*(1-t))/4.0E+0-(r*(1-s)*(1-t))/2.0E+0)/2.0E+0)-((1-s)*(1-t))/8.0E+0)+invJ.xelem(1,1)*((-((-((1-r)*(1-t2))/4.0E+0)-((1-r)*s*(1-t))/2.0E+0-((1-r2)*(1-t))/4.0E+0)/2.0E+0)-((1-r)*(1-t))/8.0E+0)+invJ.xelem(1,2)*((-((-((1-r)*(1-s)*t)/2.0E+0)-((1-r)*(1-s2))/4.0E+0-((1-r2)*(1-s))/4.0E+0)/2.0E+0)-((1-r)*(1-s))/8.0E+0);
          B.xelem(2,19) = 0;
          B.xelem(3,19) = invJ.xelem(0,0)*((-((-((1-s)*(1-t2))/4.0E+0)-((1-s2)*(1-t))/4.0E+0-(r*(1-s)*(1-t))/2.0E+0)/2.0E+0)-((1-s)*(1-t))/8.0E+0)+invJ.xelem(0,1)*((-((-((1-r)*(1-t2))/4.0E+0)-((1-r)*s*(1-t))/2.0E+0-((1-r2)*(1-t))/4.0E+0)/2.0E+0)-((1-r)*(1-t))/8.0E+0)+invJ.xelem(0,2)*((-((-((1-r)*(1-s)*t)/2.0E+0)-((1-r)*(1-s2))/4.0E+0-((1-r2)*(1-s))/4.0E+0)/2.0E+0)-((1-r)*(1-s))/8.0E+0);
          B.xelem(4,19) = invJ.xelem(2,0)*((-((-((1-s)*(1-t2))/4.0E+0)-((1-s2)*(1-t))/4.0E+0-(r*(1-s)*(1-t))/2.0E+0)/2.0E+0)-((1-s)*(1-t))/8.0E+0)+invJ.xelem(2,1)*((-((-((1-r)*(1-t2))/4.0E+0)-((1-r)*s*(1-t))/2.0E+0-((1-r2)*(1-t))/4.0E+0)/2.0E+0)-((1-r)*(1-t))/8.0E+0)+invJ.xelem(2,2)*((-((-((1-r)*(1-s)*t)/2.0E+0)-((1-r)*(1-s2))/4.0E+0-((1-r2)*(1-s))/4.0E+0)/2.0E+0)-((1-r)*(1-s))/8.0E+0);
          B.xelem(5,19) = 0;
          B.xelem(0,20) = 0;
          B.xelem(1,20) = 0;
          B.xelem(2,20) = invJ.xelem(2,0)*((-((-((1-s)*(1-t2))/4.0E+0)-((1-s2)*(1-t))/4.0E+0-(r*(1-s)*(1-t))/2.0E+0)/2.0E+0)-((1-s)*(1-t))/8.0E+0)+invJ.xelem(2,1)*((-((-((1-r)*(1-t2))/4.0E+0)-((1-r)*s*(1-t))/2.0E+0-((1-r2)*(1-t))/4.0E+0)/2.0E+0)-((1-r)*(1-t))/8.0E+0)+invJ.xelem(2,2)*((-((-((1-r)*(1-s)*t)/2.0E+0)-((1-r)*(1-s2))/4.0E+0-((1-r2)*(1-s))/4.0E+0)/2.0E+0)-((1-r)*(1-s))/8.0E+0);
          B.xelem(3,20) = 0;
          B.xelem(4,20) = invJ.xelem(1,0)*((-((-((1-s)*(1-t2))/4.0E+0)-((1-s2)*(1-t))/4.0E+0-(r*(1-s)*(1-t))/2.0E+0)/2.0E+0)-((1-s)*(1-t))/8.0E+0)+invJ.xelem(1,1)*((-((-((1-r)*(1-t2))/4.0E+0)-((1-r)*s*(1-t))/2.0E+0-((1-r2)*(1-t))/4.0E+0)/2.0E+0)-((1-r)*(1-t))/8.0E+0)+invJ.xelem(1,2)*((-((-((1-r)*(1-s)*t)/2.0E+0)-((1-r)*(1-s2))/4.0E+0-((1-r2)*(1-s))/4.0E+0)/2.0E+0)-((1-r)*(1-s))/8.0E+0);
          B.xelem(5,20) = invJ.xelem(0,0)*((-((-((1-s)*(1-t2))/4.0E+0)-((1-s2)*(1-t))/4.0E+0-(r*(1-s)*(1-t))/2.0E+0)/2.0E+0)-((1-s)*(1-t))/8.0E+0)+invJ.xelem(0,1)*((-((-((1-r)*(1-t2))/4.0E+0)-((1-r)*s*(1-t))/2.0E+0-((1-r2)*(1-t))/4.0E+0)/2.0E+0)-((1-r)*(1-t))/8.0E+0)+invJ.xelem(0,2)*((-((-((1-r)*(1-s)*t)/2.0E+0)-((1-r)*(1-s2))/4.0E+0-((1-r2)*(1-s))/4.0E+0)/2.0E+0)-((1-r)*(1-s))/8.0E+0);
          B.xelem(0,21) = invJ.xelem(0,0)*(((1-s)*(1-t))/8.0E+0-(((1-s)*(1-t2))/4.0E+0+((1-s2)*(1-t))/4.0E+0-(r*(1-s)*(1-t))/2.0E+0)/2.0E+0)+invJ.xelem(0,1)*((-((-((r+1)*(1-t2))/4.0E+0)-((r+1)*s*(1-t))/2.0E+0-((1-r2)*(1-t))/4.0E+0)/2.0E+0)-((r+1)*(1-t))/8.0E+0)+invJ.xelem(0,2)*((-((-((r+1)*(1-s)*t)/2.0E+0)-((r+1)*(1-s2))/4.0E+0-((1-r2)*(1-s))/4.0E+0)/2.0E+0)-((r+1)*(1-s))/8.0E+0);
          B.xelem(1,21) = 0;
          B.xelem(2,21) = 0;
          B.xelem(3,21) = invJ.xelem(1,0)*(((1-s)*(1-t))/8.0E+0-(((1-s)*(1-t2))/4.0E+0+((1-s2)*(1-t))/4.0E+0-(r*(1-s)*(1-t))/2.0E+0)/2.0E+0)+invJ.xelem(1,1)*((-((-((r+1)*(1-t2))/4.0E+0)-((r+1)*s*(1-t))/2.0E+0-((1-r2)*(1-t))/4.0E+0)/2.0E+0)-((r+1)*(1-t))/8.0E+0)+invJ.xelem(1,2)*((-((-((r+1)*(1-s)*t)/2.0E+0)-((r+1)*(1-s2))/4.0E+0-((1-r2)*(1-s))/4.0E+0)/2.0E+0)-((r+1)*(1-s))/8.0E+0);
          B.xelem(4,21) = 0;
          B.xelem(5,21) = invJ.xelem(2,0)*(((1-s)*(1-t))/8.0E+0-(((1-s)*(1-t2))/4.0E+0+((1-s2)*(1-t))/4.0E+0-(r*(1-s)*(1-t))/2.0E+0)/2.0E+0)+invJ.xelem(2,1)*((-((-((r+1)*(1-t2))/4.0E+0)-((r+1)*s*(1-t))/2.0E+0-((1-r2)*(1-t))/4.0E+0)/2.0E+0)-((r+1)*(1-t))/8.0E+0)+invJ.xelem(2,2)*((-((-((r+1)*(1-s)*t)/2.0E+0)-((r+1)*(1-s2))/4.0E+0-((1-r2)*(1-s))/4.0E+0)/2.0E+0)-((r+1)*(1-s))/8.0E+0);
          B.xelem(0,22) = 0;
          B.xelem(1,22) = invJ.xelem(1,0)*(((1-s)*(1-t))/8.0E+0-(((1-s)*(1-t2))/4.0E+0+((1-s2)*(1-t))/4.0E+0-(r*(1-s)*(1-t))/2.0E+0)/2.0E+0)+invJ.xelem(1,1)*((-((-((r+1)*(1-t2))/4.0E+0)-((r+1)*s*(1-t))/2.0E+0-((1-r2)*(1-t))/4.0E+0)/2.0E+0)-((r+1)*(1-t))/8.0E+0)+invJ.xelem(1,2)*((-((-((r+1)*(1-s)*t)/2.0E+0)-((r+1)*(1-s2))/4.0E+0-((1-r2)*(1-s))/4.0E+0)/2.0E+0)-((r+1)*(1-s))/8.0E+0);
          B.xelem(2,22) = 0;
          B.xelem(3,22) = invJ.xelem(0,0)*(((1-s)*(1-t))/8.0E+0-(((1-s)*(1-t2))/4.0E+0+((1-s2)*(1-t))/4.0E+0-(r*(1-s)*(1-t))/2.0E+0)/2.0E+0)+invJ.xelem(0,1)*((-((-((r+1)*(1-t2))/4.0E+0)-((r+1)*s*(1-t))/2.0E+0-((1-r2)*(1-t))/4.0E+0)/2.0E+0)-((r+1)*(1-t))/8.0E+0)+invJ.xelem(0,2)*((-((-((r+1)*(1-s)*t)/2.0E+0)-((r+1)*(1-s2))/4.0E+0-((1-r2)*(1-s))/4.0E+0)/2.0E+0)-((r+1)*(1-s))/8.0E+0);
          B.xelem(4,22) = invJ.xelem(2,0)*(((1-s)*(1-t))/8.0E+0-(((1-s)*(1-t2))/4.0E+0+((1-s2)*(1-t))/4.0E+0-(r*(1-s)*(1-t))/2.0E+0)/2.0E+0)+invJ.xelem(2,1)*((-((-((r+1)*(1-t2))/4.0E+0)-((r+1)*s*(1-t))/2.0E+0-((1-r2)*(1-t))/4.0E+0)/2.0E+0)-((r+1)*(1-t))/8.0E+0)+invJ.xelem(2,2)*((-((-((r+1)*(1-s)*t)/2.0E+0)-((r+1)*(1-s2))/4.0E+0-((1-r2)*(1-s))/4.0E+0)/2.0E+0)-((r+1)*(1-s))/8.0E+0);
          B.xelem(5,22) = 0;
          B.xelem(0,23) = 0;
          B.xelem(1,23) = 0;
          B.xelem(2,23) = invJ.xelem(2,0)*(((1-s)*(1-t))/8.0E+0-(((1-s)*(1-t2))/4.0E+0+((1-s2)*(1-t))/4.0E+0-(r*(1-s)*(1-t))/2.0E+0)/2.0E+0)+invJ.xelem(2,1)*((-((-((r+1)*(1-t2))/4.0E+0)-((r+1)*s*(1-t))/2.0E+0-((1-r2)*(1-t))/4.0E+0)/2.0E+0)-((r+1)*(1-t))/8.0E+0)+invJ.xelem(2,2)*((-((-((r+1)*(1-s)*t)/2.0E+0)-((r+1)*(1-s2))/4.0E+0-((1-r2)*(1-s))/4.0E+0)/2.0E+0)-((r+1)*(1-s))/8.0E+0);
          B.xelem(3,23) = 0;
          B.xelem(4,23) = invJ.xelem(1,0)*(((1-s)*(1-t))/8.0E+0-(((1-s)*(1-t2))/4.0E+0+((1-s2)*(1-t))/4.0E+0-(r*(1-s)*(1-t))/2.0E+0)/2.0E+0)+invJ.xelem(1,1)*((-((-((r+1)*(1-t2))/4.0E+0)-((r+1)*s*(1-t))/2.0E+0-((1-r2)*(1-t))/4.0E+0)/2.0E+0)-((r+1)*(1-t))/8.0E+0)+invJ.xelem(1,2)*((-((-((r+1)*(1-s)*t)/2.0E+0)-((r+1)*(1-s2))/4.0E+0-((1-r2)*(1-s))/4.0E+0)/2.0E+0)-((r+1)*(1-s))/8.0E+0);
          B.xelem(5,23) = invJ.xelem(0,0)*(((1-s)*(1-t))/8.0E+0-(((1-s)*(1-t2))/4.0E+0+((1-s2)*(1-t))/4.0E+0-(r*(1-s)*(1-t))/2.0E+0)/2.0E+0)+invJ.xelem(0,1)*((-((-((r+1)*(1-t2))/4.0E+0)-((r+1)*s*(1-t))/2.0E+0-((1-r2)*(1-t))/4.0E+0)/2.0E+0)-((r+1)*(1-t))/8.0E+0)+invJ.xelem(0,2)*((-((-((r+1)*(1-s)*t)/2.0E+0)-((r+1)*(1-s2))/4.0E+0-((1-r2)*(1-s))/4.0E+0)/2.0E+0)-((r+1)*(1-s))/8.0E+0);
          B.xelem(0,24) = (-(invJ.xelem(0,0)*r*(s+1)*(t+1))/2.0E+0)+(invJ.xelem(0,1)*(1-r2)*(t+1))/4.0E+0+(invJ.xelem(0,2)*(1-r2)*(s+1))/4.0E+0;
          B.xelem(1,24) = 0;
          B.xelem(2,24) = 0;
          B.xelem(3,24) = (-(invJ.xelem(1,0)*r*(s+1)*(t+1))/2.0E+0)+(invJ.xelem(1,1)*(1-r2)*(t+1))/4.0E+0+(invJ.xelem(1,2)*(1-r2)*(s+1))/4.0E+0;
          B.xelem(4,24) = 0;
          B.xelem(5,24) = (-(invJ.xelem(2,0)*r*(s+1)*(t+1))/2.0E+0)+(invJ.xelem(2,1)*(1-r2)*(t+1))/4.0E+0+(invJ.xelem(2,2)*(1-r2)*(s+1))/4.0E+0;
          B.xelem(0,25) = 0;
          B.xelem(1,25) = (-(invJ.xelem(1,0)*r*(s+1)*(t+1))/2.0E+0)+(invJ.xelem(1,1)*(1-r2)*(t+1))/4.0E+0+(invJ.xelem(1,2)*(1-r2)*(s+1))/4.0E+0;
          B.xelem(2,25) = 0;
          B.xelem(3,25) = (-(invJ.xelem(0,0)*r*(s+1)*(t+1))/2.0E+0)+(invJ.xelem(0,1)*(1-r2)*(t+1))/4.0E+0+(invJ.xelem(0,2)*(1-r2)*(s+1))/4.0E+0;
          B.xelem(4,25) = (-(invJ.xelem(2,0)*r*(s+1)*(t+1))/2.0E+0)+(invJ.xelem(2,1)*(1-r2)*(t+1))/4.0E+0+(invJ.xelem(2,2)*(1-r2)*(s+1))/4.0E+0;
          B.xelem(5,25) = 0;
          B.xelem(0,26) = 0;
          B.xelem(1,26) = 0;
          B.xelem(2,26) = (-(invJ.xelem(2,0)*r*(s+1)*(t+1))/2.0E+0)+(invJ.xelem(2,1)*(1-r2)*(t+1))/4.0E+0+(invJ.xelem(2,2)*(1-r2)*(s+1))/4.0E+0;
          B.xelem(3,26) = 0;
          B.xelem(4,26) = (-(invJ.xelem(1,0)*r*(s+1)*(t+1))/2.0E+0)+(invJ.xelem(1,1)*(1-r2)*(t+1))/4.0E+0+(invJ.xelem(1,2)*(1-r2)*(s+1))/4.0E+0;
          B.xelem(5,26) = (-(invJ.xelem(0,0)*r*(s+1)*(t+1))/2.0E+0)+(invJ.xelem(0,1)*(1-r2)*(t+1))/4.0E+0+(invJ.xelem(0,2)*(1-r2)*(s+1))/4.0E+0;
          B.xelem(0,27) = (-(invJ.xelem(0,0)*(1-s2)*(t+1))/4.0E+0)-(invJ.xelem(0,1)*(1-r)*s*(t+1))/2.0E+0+(invJ.xelem(0,2)*(1-r)*(1-s2))/4.0E+0;
          B.xelem(1,27) = 0;
          B.xelem(2,27) = 0;
          B.xelem(3,27) = (-(invJ.xelem(1,0)*(1-s2)*(t+1))/4.0E+0)-(invJ.xelem(1,1)*(1-r)*s*(t+1))/2.0E+0+(invJ.xelem(1,2)*(1-r)*(1-s2))/4.0E+0;
          B.xelem(4,27) = 0;
          B.xelem(5,27) = (-(invJ.xelem(2,0)*(1-s2)*(t+1))/4.0E+0)-(invJ.xelem(2,1)*(1-r)*s*(t+1))/2.0E+0+(invJ.xelem(2,2)*(1-r)*(1-s2))/4.0E+0;
          B.xelem(0,28) = 0;
          B.xelem(1,28) = (-(invJ.xelem(1,0)*(1-s2)*(t+1))/4.0E+0)-(invJ.xelem(1,1)*(1-r)*s*(t+1))/2.0E+0+(invJ.xelem(1,2)*(1-r)*(1-s2))/4.0E+0;
          B.xelem(2,28) = 0;
          B.xelem(3,28) = (-(invJ.xelem(0,0)*(1-s2)*(t+1))/4.0E+0)-(invJ.xelem(0,1)*(1-r)*s*(t+1))/2.0E+0+(invJ.xelem(0,2)*(1-r)*(1-s2))/4.0E+0;
          B.xelem(4,28) = (-(invJ.xelem(2,0)*(1-s2)*(t+1))/4.0E+0)-(invJ.xelem(2,1)*(1-r)*s*(t+1))/2.0E+0+(invJ.xelem(2,2)*(1-r)*(1-s2))/4.0E+0;
          B.xelem(5,28) = 0;
          B.xelem(0,29) = 0;
          B.xelem(1,29) = 0;
          B.xelem(2,29) = (-(invJ.xelem(2,0)*(1-s2)*(t+1))/4.0E+0)-(invJ.xelem(2,1)*(1-r)*s*(t+1))/2.0E+0+(invJ.xelem(2,2)*(1-r)*(1-s2))/4.0E+0;
          B.xelem(3,29) = 0;
          B.xelem(4,29) = (-(invJ.xelem(1,0)*(1-s2)*(t+1))/4.0E+0)-(invJ.xelem(1,1)*(1-r)*s*(t+1))/2.0E+0+(invJ.xelem(1,2)*(1-r)*(1-s2))/4.0E+0;
          B.xelem(5,29) = (-(invJ.xelem(0,0)*(1-s2)*(t+1))/4.0E+0)-(invJ.xelem(0,1)*(1-r)*s*(t+1))/2.0E+0+(invJ.xelem(0,2)*(1-r)*(1-s2))/4.0E+0;
          B.xelem(0,30) = (-(invJ.xelem(0,0)*r*(1-s)*(t+1))/2.0E+0)-(invJ.xelem(0,1)*(1-r2)*(t+1))/4.0E+0+(invJ.xelem(0,2)*(1-r2)*(1-s))/4.0E+0;
          B.xelem(1,30) = 0;
          B.xelem(2,30) = 0;
          B.xelem(3,30) = (-(invJ.xelem(1,0)*r*(1-s)*(t+1))/2.0E+0)-(invJ.xelem(1,1)*(1-r2)*(t+1))/4.0E+0+(invJ.xelem(1,2)*(1-r2)*(1-s))/4.0E+0;
          B.xelem(4,30) = 0;
          B.xelem(5,30) = (-(invJ.xelem(2,0)*r*(1-s)*(t+1))/2.0E+0)-(invJ.xelem(2,1)*(1-r2)*(t+1))/4.0E+0+(invJ.xelem(2,2)*(1-r2)*(1-s))/4.0E+0;
          B.xelem(0,31) = 0;
          B.xelem(1,31) = (-(invJ.xelem(1,0)*r*(1-s)*(t+1))/2.0E+0)-(invJ.xelem(1,1)*(1-r2)*(t+1))/4.0E+0+(invJ.xelem(1,2)*(1-r2)*(1-s))/4.0E+0;
          B.xelem(2,31) = 0;
          B.xelem(3,31) = (-(invJ.xelem(0,0)*r*(1-s)*(t+1))/2.0E+0)-(invJ.xelem(0,1)*(1-r2)*(t+1))/4.0E+0+(invJ.xelem(0,2)*(1-r2)*(1-s))/4.0E+0;
          B.xelem(4,31) = (-(invJ.xelem(2,0)*r*(1-s)*(t+1))/2.0E+0)-(invJ.xelem(2,1)*(1-r2)*(t+1))/4.0E+0+(invJ.xelem(2,2)*(1-r2)*(1-s))/4.0E+0;
          B.xelem(5,31) = 0;
          B.xelem(0,32) = 0;
          B.xelem(1,32) = 0;
          B.xelem(2,32) = (-(invJ.xelem(2,0)*r*(1-s)*(t+1))/2.0E+0)-(invJ.xelem(2,1)*(1-r2)*(t+1))/4.0E+0+(invJ.xelem(2,2)*(1-r2)*(1-s))/4.0E+0;
          B.xelem(3,32) = 0;
          B.xelem(4,32) = (-(invJ.xelem(1,0)*r*(1-s)*(t+1))/2.0E+0)-(invJ.xelem(1,1)*(1-r2)*(t+1))/4.0E+0+(invJ.xelem(1,2)*(1-r2)*(1-s))/4.0E+0;
          B.xelem(5,32) = (-(invJ.xelem(0,0)*r*(1-s)*(t+1))/2.0E+0)-(invJ.xelem(0,1)*(1-r2)*(t+1))/4.0E+0+(invJ.xelem(0,2)*(1-r2)*(1-s))/4.0E+0;
          B.xelem(0,33) = (invJ.xelem(0,0)*(1-s2)*(t+1))/4.0E+0-(invJ.xelem(0,1)*(r+1)*s*(t+1))/2.0E+0+(invJ.xelem(0,2)*(r+1)*(1-s2))/4.0E+0;
          B.xelem(1,33) = 0;
          B.xelem(2,33) = 0;
          B.xelem(3,33) = (invJ.xelem(1,0)*(1-s2)*(t+1))/4.0E+0-(invJ.xelem(1,1)*(r+1)*s*(t+1))/2.0E+0+(invJ.xelem(1,2)*(r+1)*(1-s2))/4.0E+0;
          B.xelem(4,33) = 0;
          B.xelem(5,33) = (invJ.xelem(2,0)*(1-s2)*(t+1))/4.0E+0-(invJ.xelem(2,1)*(r+1)*s*(t+1))/2.0E+0+(invJ.xelem(2,2)*(r+1)*(1-s2))/4.0E+0;
          B.xelem(0,34) = 0;
          B.xelem(1,34) = (invJ.xelem(1,0)*(1-s2)*(t+1))/4.0E+0-(invJ.xelem(1,1)*(r+1)*s*(t+1))/2.0E+0+(invJ.xelem(1,2)*(r+1)*(1-s2))/4.0E+0;
          B.xelem(2,34) = 0;
          B.xelem(3,34) = (invJ.xelem(0,0)*(1-s2)*(t+1))/4.0E+0-(invJ.xelem(0,1)*(r+1)*s*(t+1))/2.0E+0+(invJ.xelem(0,2)*(r+1)*(1-s2))/4.0E+0;
          B.xelem(4,34) = (invJ.xelem(2,0)*(1-s2)*(t+1))/4.0E+0-(invJ.xelem(2,1)*(r+1)*s*(t+1))/2.0E+0+(invJ.xelem(2,2)*(r+1)*(1-s2))/4.0E+0;
          B.xelem(5,34) = 0;
          B.xelem(0,35) = 0;
          B.xelem(1,35) = 0;
          B.xelem(2,35) = (invJ.xelem(2,0)*(1-s2)*(t+1))/4.0E+0-(invJ.xelem(2,1)*(r+1)*s*(t+1))/2.0E+0+(invJ.xelem(2,2)*(r+1)*(1-s2))/4.0E+0;
          B.xelem(3,35) = 0;
          B.xelem(4,35) = (invJ.xelem(1,0)*(1-s2)*(t+1))/4.0E+0-(invJ.xelem(1,1)*(r+1)*s*(t+1))/2.0E+0+(invJ.xelem(1,2)*(r+1)*(1-s2))/4.0E+0;
          B.xelem(5,35) = (invJ.xelem(0,0)*(1-s2)*(t+1))/4.0E+0-(invJ.xelem(0,1)*(r+1)*s*(t+1))/2.0E+0+(invJ.xelem(0,2)*(r+1)*(1-s2))/4.0E+0;
          B.xelem(0,36) = (-(invJ.xelem(0,0)*r*(s+1)*(1-t))/2.0E+0)+(invJ.xelem(0,1)*(1-r2)*(1-t))/4.0E+0-(invJ.xelem(0,2)*(1-r2)*(s+1))/4.0E+0;
          B.xelem(1,36) = 0;
          B.xelem(2,36) = 0;
          B.xelem(3,36) = (-(invJ.xelem(1,0)*r*(s+1)*(1-t))/2.0E+0)+(invJ.xelem(1,1)*(1-r2)*(1-t))/4.0E+0-(invJ.xelem(1,2)*(1-r2)*(s+1))/4.0E+0;
          B.xelem(4,36) = 0;
          B.xelem(5,36) = (-(invJ.xelem(2,0)*r*(s+1)*(1-t))/2.0E+0)+(invJ.xelem(2,1)*(1-r2)*(1-t))/4.0E+0-(invJ.xelem(2,2)*(1-r2)*(s+1))/4.0E+0;
          B.xelem(0,37) = 0;
          B.xelem(1,37) = (-(invJ.xelem(1,0)*r*(s+1)*(1-t))/2.0E+0)+(invJ.xelem(1,1)*(1-r2)*(1-t))/4.0E+0-(invJ.xelem(1,2)*(1-r2)*(s+1))/4.0E+0;
          B.xelem(2,37) = 0;
          B.xelem(3,37) = (-(invJ.xelem(0,0)*r*(s+1)*(1-t))/2.0E+0)+(invJ.xelem(0,1)*(1-r2)*(1-t))/4.0E+0-(invJ.xelem(0,2)*(1-r2)*(s+1))/4.0E+0;
          B.xelem(4,37) = (-(invJ.xelem(2,0)*r*(s+1)*(1-t))/2.0E+0)+(invJ.xelem(2,1)*(1-r2)*(1-t))/4.0E+0-(invJ.xelem(2,2)*(1-r2)*(s+1))/4.0E+0;
          B.xelem(5,37) = 0;
          B.xelem(0,38) = 0;
          B.xelem(1,38) = 0;
          B.xelem(2,38) = (-(invJ.xelem(2,0)*r*(s+1)*(1-t))/2.0E+0)+(invJ.xelem(2,1)*(1-r2)*(1-t))/4.0E+0-(invJ.xelem(2,2)*(1-r2)*(s+1))/4.0E+0;
          B.xelem(3,38) = 0;
          B.xelem(4,38) = (-(invJ.xelem(1,0)*r*(s+1)*(1-t))/2.0E+0)+(invJ.xelem(1,1)*(1-r2)*(1-t))/4.0E+0-(invJ.xelem(1,2)*(1-r2)*(s+1))/4.0E+0;
          B.xelem(5,38) = (-(invJ.xelem(0,0)*r*(s+1)*(1-t))/2.0E+0)+(invJ.xelem(0,1)*(1-r2)*(1-t))/4.0E+0-(invJ.xelem(0,2)*(1-r2)*(s+1))/4.0E+0;
          B.xelem(0,39) = (-(invJ.xelem(0,0)*(1-s2)*(1-t))/4.0E+0)-(invJ.xelem(0,1)*(1-r)*s*(1-t))/2.0E+0-(invJ.xelem(0,2)*(1-r)*(1-s2))/4.0E+0;
          B.xelem(1,39) = 0;
          B.xelem(2,39) = 0;
          B.xelem(3,39) = (-(invJ.xelem(1,0)*(1-s2)*(1-t))/4.0E+0)-(invJ.xelem(1,1)*(1-r)*s*(1-t))/2.0E+0-(invJ.xelem(1,2)*(1-r)*(1-s2))/4.0E+0;
          B.xelem(4,39) = 0;
          B.xelem(5,39) = (-(invJ.xelem(2,0)*(1-s2)*(1-t))/4.0E+0)-(invJ.xelem(2,1)*(1-r)*s*(1-t))/2.0E+0-(invJ.xelem(2,2)*(1-r)*(1-s2))/4.0E+0;
          B.xelem(0,40) = 0;
          B.xelem(1,40) = (-(invJ.xelem(1,0)*(1-s2)*(1-t))/4.0E+0)-(invJ.xelem(1,1)*(1-r)*s*(1-t))/2.0E+0-(invJ.xelem(1,2)*(1-r)*(1-s2))/4.0E+0;
          B.xelem(2,40) = 0;
          B.xelem(3,40) = (-(invJ.xelem(0,0)*(1-s2)*(1-t))/4.0E+0)-(invJ.xelem(0,1)*(1-r)*s*(1-t))/2.0E+0-(invJ.xelem(0,2)*(1-r)*(1-s2))/4.0E+0;
          B.xelem(4,40) = (-(invJ.xelem(2,0)*(1-s2)*(1-t))/4.0E+0)-(invJ.xelem(2,1)*(1-r)*s*(1-t))/2.0E+0-(invJ.xelem(2,2)*(1-r)*(1-s2))/4.0E+0;
          B.xelem(5,40) = 0;
          B.xelem(0,41) = 0;
          B.xelem(1,41) = 0;
          B.xelem(2,41) = (-(invJ.xelem(2,0)*(1-s2)*(1-t))/4.0E+0)-(invJ.xelem(2,1)*(1-r)*s*(1-t))/2.0E+0-(invJ.xelem(2,2)*(1-r)*(1-s2))/4.0E+0;
          B.xelem(3,41) = 0;
          B.xelem(4,41) = (-(invJ.xelem(1,0)*(1-s2)*(1-t))/4.0E+0)-(invJ.xelem(1,1)*(1-r)*s*(1-t))/2.0E+0-(invJ.xelem(1,2)*(1-r)*(1-s2))/4.0E+0;
          B.xelem(5,41) = (-(invJ.xelem(0,0)*(1-s2)*(1-t))/4.0E+0)-(invJ.xelem(0,1)*(1-r)*s*(1-t))/2.0E+0-(invJ.xelem(0,2)*(1-r)*(1-s2))/4.0E+0;
          B.xelem(0,42) = (-(invJ.xelem(0,0)*r*(1-s)*(1-t))/2.0E+0)-(invJ.xelem(0,1)*(1-r2)*(1-t))/4.0E+0-(invJ.xelem(0,2)*(1-r2)*(1-s))/4.0E+0;
          B.xelem(1,42) = 0;
          B.xelem(2,42) = 0;
          B.xelem(3,42) = (-(invJ.xelem(1,0)*r*(1-s)*(1-t))/2.0E+0)-(invJ.xelem(1,1)*(1-r2)*(1-t))/4.0E+0-(invJ.xelem(1,2)*(1-r2)*(1-s))/4.0E+0;
          B.xelem(4,42) = 0;
          B.xelem(5,42) = (-(invJ.xelem(2,0)*r*(1-s)*(1-t))/2.0E+0)-(invJ.xelem(2,1)*(1-r2)*(1-t))/4.0E+0-(invJ.xelem(2,2)*(1-r2)*(1-s))/4.0E+0;
          B.xelem(0,43) = 0;
          B.xelem(1,43) = (-(invJ.xelem(1,0)*r*(1-s)*(1-t))/2.0E+0)-(invJ.xelem(1,1)*(1-r2)*(1-t))/4.0E+0-(invJ.xelem(1,2)*(1-r2)*(1-s))/4.0E+0;
          B.xelem(2,43) = 0;
          B.xelem(3,43) = (-(invJ.xelem(0,0)*r*(1-s)*(1-t))/2.0E+0)-(invJ.xelem(0,1)*(1-r2)*(1-t))/4.0E+0-(invJ.xelem(0,2)*(1-r2)*(1-s))/4.0E+0;
          B.xelem(4,43) = (-(invJ.xelem(2,0)*r*(1-s)*(1-t))/2.0E+0)-(invJ.xelem(2,1)*(1-r2)*(1-t))/4.0E+0-(invJ.xelem(2,2)*(1-r2)*(1-s))/4.0E+0;
          B.xelem(5,43) = 0;
          B.xelem(0,44) = 0;
          B.xelem(1,44) = 0;
          B.xelem(2,44) = (-(invJ.xelem(2,0)*r*(1-s)*(1-t))/2.0E+0)-(invJ.xelem(2,1)*(1-r2)*(1-t))/4.0E+0-(invJ.xelem(2,2)*(1-r2)*(1-s))/4.0E+0;
          B.xelem(3,44) = 0;
          B.xelem(4,44) = (-(invJ.xelem(1,0)*r*(1-s)*(1-t))/2.0E+0)-(invJ.xelem(1,1)*(1-r2)*(1-t))/4.0E+0-(invJ.xelem(1,2)*(1-r2)*(1-s))/4.0E+0;
          B.xelem(5,44) = (-(invJ.xelem(0,0)*r*(1-s)*(1-t))/2.0E+0)-(invJ.xelem(0,1)*(1-r2)*(1-t))/4.0E+0-(invJ.xelem(0,2)*(1-r2)*(1-s))/4.0E+0;
          B.xelem(0,45) = (invJ.xelem(0,0)*(1-s2)*(1-t))/4.0E+0-(invJ.xelem(0,1)*(r+1)*s*(1-t))/2.0E+0-(invJ.xelem(0,2)*(r+1)*(1-s2))/4.0E+0;
          B.xelem(1,45) = 0;
          B.xelem(2,45) = 0;
          B.xelem(3,45) = (invJ.xelem(1,0)*(1-s2)*(1-t))/4.0E+0-(invJ.xelem(1,1)*(r+1)*s*(1-t))/2.0E+0-(invJ.xelem(1,2)*(r+1)*(1-s2))/4.0E+0;
          B.xelem(4,45) = 0;
          B.xelem(5,45) = (invJ.xelem(2,0)*(1-s2)*(1-t))/4.0E+0-(invJ.xelem(2,1)*(r+1)*s*(1-t))/2.0E+0-(invJ.xelem(2,2)*(r+1)*(1-s2))/4.0E+0;
          B.xelem(0,46) = 0;
          B.xelem(1,46) = (invJ.xelem(1,0)*(1-s2)*(1-t))/4.0E+0-(invJ.xelem(1,1)*(r+1)*s*(1-t))/2.0E+0-(invJ.xelem(1,2)*(r+1)*(1-s2))/4.0E+0;
          B.xelem(2,46) = 0;
          B.xelem(3,46) = (invJ.xelem(0,0)*(1-s2)*(1-t))/4.0E+0-(invJ.xelem(0,1)*(r+1)*s*(1-t))/2.0E+0-(invJ.xelem(0,2)*(r+1)*(1-s2))/4.0E+0;
          B.xelem(4,46) = (invJ.xelem(2,0)*(1-s2)*(1-t))/4.0E+0-(invJ.xelem(2,1)*(r+1)*s*(1-t))/2.0E+0-(invJ.xelem(2,2)*(r+1)*(1-s2))/4.0E+0;
          B.xelem(5,46) = 0;
          B.xelem(0,47) = 0;
          B.xelem(1,47) = 0;
          B.xelem(2,47) = (invJ.xelem(2,0)*(1-s2)*(1-t))/4.0E+0-(invJ.xelem(2,1)*(r+1)*s*(1-t))/2.0E+0-(invJ.xelem(2,2)*(r+1)*(1-s2))/4.0E+0;
          B.xelem(3,47) = 0;
          B.xelem(4,47) = (invJ.xelem(1,0)*(1-s2)*(1-t))/4.0E+0-(invJ.xelem(1,1)*(r+1)*s*(1-t))/2.0E+0-(invJ.xelem(1,2)*(r+1)*(1-s2))/4.0E+0;
          B.xelem(5,47) = (invJ.xelem(0,0)*(1-s2)*(1-t))/4.0E+0-(invJ.xelem(0,1)*(r+1)*s*(1-t))/2.0E+0-(invJ.xelem(0,2)*(r+1)*(1-s2))/4.0E+0;
          B.xelem(0,48) = (invJ.xelem(0,0)*(s+1)*(1-t2))/4.0E+0+(invJ.xelem(0,1)*(r+1)*(1-t2))/4.0E+0-(invJ.xelem(0,2)*(r+1)*(s+1)*t)/2.0E+0;
          B.xelem(1,48) = 0;
          B.xelem(2,48) = 0;
          B.xelem(3,48) = (invJ.xelem(1,0)*(s+1)*(1-t2))/4.0E+0+(invJ.xelem(1,1)*(r+1)*(1-t2))/4.0E+0-(invJ.xelem(1,2)*(r+1)*(s+1)*t)/2.0E+0;
          B.xelem(4,48) = 0;
          B.xelem(5,48) = (invJ.xelem(2,0)*(s+1)*(1-t2))/4.0E+0+(invJ.xelem(2,1)*(r+1)*(1-t2))/4.0E+0-(invJ.xelem(2,2)*(r+1)*(s+1)*t)/2.0E+0;
          B.xelem(0,49) = 0;
          B.xelem(1,49) = (invJ.xelem(1,0)*(s+1)*(1-t2))/4.0E+0+(invJ.xelem(1,1)*(r+1)*(1-t2))/4.0E+0-(invJ.xelem(1,2)*(r+1)*(s+1)*t)/2.0E+0;
          B.xelem(2,49) = 0;
          B.xelem(3,49) = (invJ.xelem(0,0)*(s+1)*(1-t2))/4.0E+0+(invJ.xelem(0,1)*(r+1)*(1-t2))/4.0E+0-(invJ.xelem(0,2)*(r+1)*(s+1)*t)/2.0E+0;
          B.xelem(4,49) = (invJ.xelem(2,0)*(s+1)*(1-t2))/4.0E+0+(invJ.xelem(2,1)*(r+1)*(1-t2))/4.0E+0-(invJ.xelem(2,2)*(r+1)*(s+1)*t)/2.0E+0;
          B.xelem(5,49) = 0;
          B.xelem(0,50) = 0;
          B.xelem(1,50) = 0;
          B.xelem(2,50) = (invJ.xelem(2,0)*(s+1)*(1-t2))/4.0E+0+(invJ.xelem(2,1)*(r+1)*(1-t2))/4.0E+0-(invJ.xelem(2,2)*(r+1)*(s+1)*t)/2.0E+0;
          B.xelem(3,50) = 0;
          B.xelem(4,50) = (invJ.xelem(1,0)*(s+1)*(1-t2))/4.0E+0+(invJ.xelem(1,1)*(r+1)*(1-t2))/4.0E+0-(invJ.xelem(1,2)*(r+1)*(s+1)*t)/2.0E+0;
          B.xelem(5,50) = (invJ.xelem(0,0)*(s+1)*(1-t2))/4.0E+0+(invJ.xelem(0,1)*(r+1)*(1-t2))/4.0E+0-(invJ.xelem(0,2)*(r+1)*(s+1)*t)/2.0E+0;
          B.xelem(0,51) = (-(invJ.xelem(0,0)*(s+1)*(1-t2))/4.0E+0)+(invJ.xelem(0,1)*(1-r)*(1-t2))/4.0E+0-(invJ.xelem(0,2)*(1-r)*(s+1)*t)/2.0E+0;
          B.xelem(1,51) = 0;
          B.xelem(2,51) = 0;
          B.xelem(3,51) = (-(invJ.xelem(1,0)*(s+1)*(1-t2))/4.0E+0)+(invJ.xelem(1,1)*(1-r)*(1-t2))/4.0E+0-(invJ.xelem(1,2)*(1-r)*(s+1)*t)/2.0E+0;
          B.xelem(4,51) = 0;
          B.xelem(5,51) = (-(invJ.xelem(2,0)*(s+1)*(1-t2))/4.0E+0)+(invJ.xelem(2,1)*(1-r)*(1-t2))/4.0E+0-(invJ.xelem(2,2)*(1-r)*(s+1)*t)/2.0E+0;
          B.xelem(0,52) = 0;
          B.xelem(1,52) = (-(invJ.xelem(1,0)*(s+1)*(1-t2))/4.0E+0)+(invJ.xelem(1,1)*(1-r)*(1-t2))/4.0E+0-(invJ.xelem(1,2)*(1-r)*(s+1)*t)/2.0E+0;
          B.xelem(2,52) = 0;
          B.xelem(3,52) = (-(invJ.xelem(0,0)*(s+1)*(1-t2))/4.0E+0)+(invJ.xelem(0,1)*(1-r)*(1-t2))/4.0E+0-(invJ.xelem(0,2)*(1-r)*(s+1)*t)/2.0E+0;
          B.xelem(4,52) = (-(invJ.xelem(2,0)*(s+1)*(1-t2))/4.0E+0)+(invJ.xelem(2,1)*(1-r)*(1-t2))/4.0E+0-(invJ.xelem(2,2)*(1-r)*(s+1)*t)/2.0E+0;
          B.xelem(5,52) = 0;
          B.xelem(0,53) = 0;
          B.xelem(1,53) = 0;
          B.xelem(2,53) = (-(invJ.xelem(2,0)*(s+1)*(1-t2))/4.0E+0)+(invJ.xelem(2,1)*(1-r)*(1-t2))/4.0E+0-(invJ.xelem(2,2)*(1-r)*(s+1)*t)/2.0E+0;
          B.xelem(3,53) = 0;
          B.xelem(4,53) = (-(invJ.xelem(1,0)*(s+1)*(1-t2))/4.0E+0)+(invJ.xelem(1,1)*(1-r)*(1-t2))/4.0E+0-(invJ.xelem(1,2)*(1-r)*(s+1)*t)/2.0E+0;
          B.xelem(5,53) = (-(invJ.xelem(0,0)*(s+1)*(1-t2))/4.0E+0)+(invJ.xelem(0,1)*(1-r)*(1-t2))/4.0E+0-(invJ.xelem(0,2)*(1-r)*(s+1)*t)/2.0E+0;
          B.xelem(0,54) = (-(invJ.xelem(0,0)*(1-s)*(1-t2))/4.0E+0)-(invJ.xelem(0,1)*(1-r)*(1-t2))/4.0E+0-(invJ.xelem(0,2)*(1-r)*(1-s)*t)/2.0E+0;
          B.xelem(1,54) = 0;
          B.xelem(2,54) = 0;
          B.xelem(3,54) = (-(invJ.xelem(1,0)*(1-s)*(1-t2))/4.0E+0)-(invJ.xelem(1,1)*(1-r)*(1-t2))/4.0E+0-(invJ.xelem(1,2)*(1-r)*(1-s)*t)/2.0E+0;
          B.xelem(4,54) = 0;
          B.xelem(5,54) = (-(invJ.xelem(2,0)*(1-s)*(1-t2))/4.0E+0)-(invJ.xelem(2,1)*(1-r)*(1-t2))/4.0E+0-(invJ.xelem(2,2)*(1-r)*(1-s)*t)/2.0E+0;
          B.xelem(0,55) = 0;
          B.xelem(1,55) = (-(invJ.xelem(1,0)*(1-s)*(1-t2))/4.0E+0)-(invJ.xelem(1,1)*(1-r)*(1-t2))/4.0E+0-(invJ.xelem(1,2)*(1-r)*(1-s)*t)/2.0E+0;
          B.xelem(2,55) = 0;
          B.xelem(3,55) = (-(invJ.xelem(0,0)*(1-s)*(1-t2))/4.0E+0)-(invJ.xelem(0,1)*(1-r)*(1-t2))/4.0E+0-(invJ.xelem(0,2)*(1-r)*(1-s)*t)/2.0E+0;
          B.xelem(4,55) = (-(invJ.xelem(2,0)*(1-s)*(1-t2))/4.0E+0)-(invJ.xelem(2,1)*(1-r)*(1-t2))/4.0E+0-(invJ.xelem(2,2)*(1-r)*(1-s)*t)/2.0E+0;
          B.xelem(5,55) = 0;
          B.xelem(0,56) = 0;
          B.xelem(1,56) = 0;
          B.xelem(2,56) = (-(invJ.xelem(2,0)*(1-s)*(1-t2))/4.0E+0)-(invJ.xelem(2,1)*(1-r)*(1-t2))/4.0E+0-(invJ.xelem(2,2)*(1-r)*(1-s)*t)/2.0E+0;
          B.xelem(3,56) = 0;
          B.xelem(4,56) = (-(invJ.xelem(1,0)*(1-s)*(1-t2))/4.0E+0)-(invJ.xelem(1,1)*(1-r)*(1-t2))/4.0E+0-(invJ.xelem(1,2)*(1-r)*(1-s)*t)/2.0E+0;
          B.xelem(5,56) = (-(invJ.xelem(0,0)*(1-s)*(1-t2))/4.0E+0)-(invJ.xelem(0,1)*(1-r)*(1-t2))/4.0E+0-(invJ.xelem(0,2)*(1-r)*(1-s)*t)/2.0E+0;
          B.xelem(0,57) = (invJ.xelem(0,0)*(1-s)*(1-t2))/4.0E+0-(invJ.xelem(0,1)*(r+1)*(1-t2))/4.0E+0-(invJ.xelem(0,2)*(r+1)*(1-s)*t)/2.0E+0;
          B.xelem(1,57) = 0;
          B.xelem(2,57) = 0;
          B.xelem(3,57) = (invJ.xelem(1,0)*(1-s)*(1-t2))/4.0E+0-(invJ.xelem(1,1)*(r+1)*(1-t2))/4.0E+0-(invJ.xelem(1,2)*(r+1)*(1-s)*t)/2.0E+0;
          B.xelem(4,57) = 0;
          B.xelem(5,57) = (invJ.xelem(2,0)*(1-s)*(1-t2))/4.0E+0-(invJ.xelem(2,1)*(r+1)*(1-t2))/4.0E+0-(invJ.xelem(2,2)*(r+1)*(1-s)*t)/2.0E+0;
          B.xelem(0,58) = 0;
          B.xelem(1,58) = (invJ.xelem(1,0)*(1-s)*(1-t2))/4.0E+0-(invJ.xelem(1,1)*(r+1)*(1-t2))/4.0E+0-(invJ.xelem(1,2)*(r+1)*(1-s)*t)/2.0E+0;
          B.xelem(2,58) = 0;
          B.xelem(3,58) = (invJ.xelem(0,0)*(1-s)*(1-t2))/4.0E+0-(invJ.xelem(0,1)*(r+1)*(1-t2))/4.0E+0-(invJ.xelem(0,2)*(r+1)*(1-s)*t)/2.0E+0;
          B.xelem(4,58) = (invJ.xelem(2,0)*(1-s)*(1-t2))/4.0E+0-(invJ.xelem(2,1)*(r+1)*(1-t2))/4.0E+0-(invJ.xelem(2,2)*(r+1)*(1-s)*t)/2.0E+0;
          B.xelem(5,58) = 0;
          B.xelem(0,59) = 0;
          B.xelem(1,59) = 0;
          B.xelem(2,59) = (invJ.xelem(2,0)*(1-s)*(1-t2))/4.0E+0-(invJ.xelem(2,1)*(r+1)*(1-t2))/4.0E+0-(invJ.xelem(2,2)*(r+1)*(1-s)*t)/2.0E+0;
          B.xelem(3,59) = 0;
          B.xelem(4,59) = (invJ.xelem(1,0)*(1-s)*(1-t2))/4.0E+0-(invJ.xelem(1,1)*(r+1)*(1-t2))/4.0E+0-(invJ.xelem(1,2)*(r+1)*(1-s)*t)/2.0E+0;
          B.xelem(5,59) = (invJ.xelem(0,0)*(1-s)*(1-t2))/4.0E+0-(invJ.xelem(0,1)*(r+1)*(1-t2))/4.0E+0-(invJ.xelem(0,2)*(r+1)*(1-s)*t)/2.0E+0;
     }

     virtual Matrix InterpGaussToNodal(FemMatrixType eMatType, const Matrix& taug) const {
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

private:
     void ScalarInterpMatrix(const ColumnVector& rv, Matrix& Hs, octave_idx_type irow) const {
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

          Hs.xelem(irow, 0) = ((r+1)*(s+1)*(t+1))/8.0E+0-(((r+1)*(s+1)*(1-t2))/4.0E+0+((r+1)*(1-s2)*(t+1))/4.0E+0+((1-r2)*(s+1)*(t+1))/4.0E+0)/2.0E+0;
          Hs.xelem(irow, 1) = ((1-r)*(s+1)*(t+1))/8.0E+0-(((1-r)*(s+1)*(1-t2))/4.0E+0+((1-r)*(1-s2)*(t+1))/4.0E+0+((1-r2)*(s+1)*(t+1))/4.0E+0)/2.0E+0;
          Hs.xelem(irow, 2) = ((1-r)*(1-s)*(t+1))/8.0E+0-(((1-r)*(1-s)*(1-t2))/4.0E+0+((1-r)*(1-s2)*(t+1))/4.0E+0+((1-r2)*(1-s)*(t+1))/4.0E+0)/2.0E+0;
          Hs.xelem(irow, 3) = ((r+1)*(1-s)*(t+1))/8.0E+0-(((r+1)*(1-s)*(1-t2))/4.0E+0+((r+1)*(1-s2)*(t+1))/4.0E+0+((1-r2)*(1-s)*(t+1))/4.0E+0)/2.0E+0;
          Hs.xelem(irow, 4) = ((r+1)*(s+1)*(1-t))/8.0E+0-(((r+1)*(s+1)*(1-t2))/4.0E+0+((r+1)*(1-s2)*(1-t))/4.0E+0+((1-r2)*(s+1)*(1-t))/4.0E+0)/2.0E+0;
          Hs.xelem(irow, 5) = ((1-r)*(s+1)*(1-t))/8.0E+0-(((1-r)*(s+1)*(1-t2))/4.0E+0+((1-r)*(1-s2)*(1-t))/4.0E+0+((1-r2)*(s+1)*(1-t))/4.0E+0)/2.0E+0;
          Hs.xelem(irow, 6) = ((1-r)*(1-s)*(1-t))/8.0E+0-(((1-r)*(1-s)*(1-t2))/4.0E+0+((1-r)*(1-s2)*(1-t))/4.0E+0+((1-r2)*(1-s)*(1-t))/4.0E+0)/2.0E+0;
          Hs.xelem(irow, 7) = ((r+1)*(1-s)*(1-t))/8.0E+0-(((r+1)*(1-s)*(1-t2))/4.0E+0+((r+1)*(1-s2)*(1-t))/4.0E+0+((1-r2)*(1-s)*(1-t))/4.0E+0)/2.0E+0;
          Hs.xelem(irow, 8) = ((1-r2)*(s+1)*(t+1))/4.0E+0;
          Hs.xelem(irow, 9) = ((1-r)*(1-s2)*(t+1))/4.0E+0;
          Hs.xelem(irow, 10) = ((1-r2)*(1-s)*(t+1))/4.0E+0;
          Hs.xelem(irow, 11) = ((r+1)*(1-s2)*(t+1))/4.0E+0;
          Hs.xelem(irow, 12) = ((1-r2)*(s+1)*(1-t))/4.0E+0;
          Hs.xelem(irow, 13) = ((1-r)*(1-s2)*(1-t))/4.0E+0;
          Hs.xelem(irow, 14) = ((1-r2)*(1-s)*(1-t))/4.0E+0;
          Hs.xelem(irow, 15) = ((r+1)*(1-s2)*(1-t))/4.0E+0;
          Hs.xelem(irow, 16) = ((r+1)*(s+1)*(1-t2))/4.0E+0;
          Hs.xelem(irow, 17) = ((1-r)*(s+1)*(1-t2))/4.0E+0;
          Hs.xelem(irow, 18) = ((1-r)*(1-s)*(1-t2))/4.0E+0;
          Hs.xelem(irow, 19) = ((r+1)*(1-s)*(1-t2))/4.0E+0;
     }

     virtual void TemperatureGradientMatrix(const ColumnVector& rv, const Matrix& invJ, Matrix& Bt) const {
          FEM_ASSERT(Bt.rows() == 3);
          FEM_ASSERT(Bt.columns() == 20);
          
          const double r = rv.xelem(0);
          const double s = rv.xelem(1);
          const double t = rv.xelem(2);
          const double r_2 = r * r;
          const double s_2 = s * s;
          const double t_2 = t * t;
          
          Bt.xelem(0,0) = invJ.xelem(0,0)*(((s+1)*(t+1))/8.0E+0-(((s+1)*(1-t_2))/4.0E+0+((1-s_2)*(t+1))/4.0E+0-(r*(s+1)*(t+1))/2.0E+0)/2.0E+0)+invJ.xelem(0,1)*(((r+1)*(t+1))/8.0E+0-(((r+1)*(1-t_2))/4.0E+0-((r+1)*s*(t+1))/2.0E+0+((1-r_2)*(t+1))/4.0E+0)/2.0E+0)+invJ.xelem(0,2)*(((r+1)*(s+1))/8.0E+0-((-((r+1)*(s+1)*t)/2.0E+0)+((r+1)*(1-s_2))/4.0E+0+((1-r_2)*(s+1))/4.0E+0)/2.0E+0);
          Bt.xelem(1,0) = invJ.xelem(1,0)*(((s+1)*(t+1))/8.0E+0-(((s+1)*(1-t_2))/4.0E+0+((1-s_2)*(t+1))/4.0E+0-(r*(s+1)*(t+1))/2.0E+0)/2.0E+0)+invJ.xelem(1,1)*(((r+1)*(t+1))/8.0E+0-(((r+1)*(1-t_2))/4.0E+0-((r+1)*s*(t+1))/2.0E+0+((1-r_2)*(t+1))/4.0E+0)/2.0E+0)+invJ.xelem(1,2)*(((r+1)*(s+1))/8.0E+0-((-((r+1)*(s+1)*t)/2.0E+0)+((r+1)*(1-s_2))/4.0E+0+((1-r_2)*(s+1))/4.0E+0)/2.0E+0);
          Bt.xelem(2,0) = invJ.xelem(2,0)*(((s+1)*(t+1))/8.0E+0-(((s+1)*(1-t_2))/4.0E+0+((1-s_2)*(t+1))/4.0E+0-(r*(s+1)*(t+1))/2.0E+0)/2.0E+0)+invJ.xelem(2,1)*(((r+1)*(t+1))/8.0E+0-(((r+1)*(1-t_2))/4.0E+0-((r+1)*s*(t+1))/2.0E+0+((1-r_2)*(t+1))/4.0E+0)/2.0E+0)+invJ.xelem(2,2)*(((r+1)*(s+1))/8.0E+0-((-((r+1)*(s+1)*t)/2.0E+0)+((r+1)*(1-s_2))/4.0E+0+((1-r_2)*(s+1))/4.0E+0)/2.0E+0);
          Bt.xelem(0,1) = invJ.xelem(0,0)*((-((-((s+1)*(1-t_2))/4.0E+0)-((1-s_2)*(t+1))/4.0E+0-(r*(s+1)*(t+1))/2.0E+0)/2.0E+0)-((s+1)*(t+1))/8.0E+0)+invJ.xelem(0,1)*(((1-r)*(t+1))/8.0E+0-(((1-r)*(1-t_2))/4.0E+0-((1-r)*s*(t+1))/2.0E+0+((1-r_2)*(t+1))/4.0E+0)/2.0E+0)+invJ.xelem(0,2)*(((1-r)*(s+1))/8.0E+0-((-((1-r)*(s+1)*t)/2.0E+0)+((1-r)*(1-s_2))/4.0E+0+((1-r_2)*(s+1))/4.0E+0)/2.0E+0);
          Bt.xelem(1,1) = invJ.xelem(1,0)*((-((-((s+1)*(1-t_2))/4.0E+0)-((1-s_2)*(t+1))/4.0E+0-(r*(s+1)*(t+1))/2.0E+0)/2.0E+0)-((s+1)*(t+1))/8.0E+0)+invJ.xelem(1,1)*(((1-r)*(t+1))/8.0E+0-(((1-r)*(1-t_2))/4.0E+0-((1-r)*s*(t+1))/2.0E+0+((1-r_2)*(t+1))/4.0E+0)/2.0E+0)+invJ.xelem(1,2)*(((1-r)*(s+1))/8.0E+0-((-((1-r)*(s+1)*t)/2.0E+0)+((1-r)*(1-s_2))/4.0E+0+((1-r_2)*(s+1))/4.0E+0)/2.0E+0);
          Bt.xelem(2,1) = invJ.xelem(2,0)*((-((-((s+1)*(1-t_2))/4.0E+0)-((1-s_2)*(t+1))/4.0E+0-(r*(s+1)*(t+1))/2.0E+0)/2.0E+0)-((s+1)*(t+1))/8.0E+0)+invJ.xelem(2,1)*(((1-r)*(t+1))/8.0E+0-(((1-r)*(1-t_2))/4.0E+0-((1-r)*s*(t+1))/2.0E+0+((1-r_2)*(t+1))/4.0E+0)/2.0E+0)+invJ.xelem(2,2)*(((1-r)*(s+1))/8.0E+0-((-((1-r)*(s+1)*t)/2.0E+0)+((1-r)*(1-s_2))/4.0E+0+((1-r_2)*(s+1))/4.0E+0)/2.0E+0);
          Bt.xelem(0,2) = invJ.xelem(0,0)*((-((-((1-s)*(1-t_2))/4.0E+0)-((1-s_2)*(t+1))/4.0E+0-(r*(1-s)*(t+1))/2.0E+0)/2.0E+0)-((1-s)*(t+1))/8.0E+0)+invJ.xelem(0,1)*((-((-((1-r)*(1-t_2))/4.0E+0)-((1-r)*s*(t+1))/2.0E+0-((1-r_2)*(t+1))/4.0E+0)/2.0E+0)-((1-r)*(t+1))/8.0E+0)+invJ.xelem(0,2)*(((1-r)*(1-s))/8.0E+0-((-((1-r)*(1-s)*t)/2.0E+0)+((1-r)*(1-s_2))/4.0E+0+((1-r_2)*(1-s))/4.0E+0)/2.0E+0);
          Bt.xelem(1,2) = invJ.xelem(1,0)*((-((-((1-s)*(1-t_2))/4.0E+0)-((1-s_2)*(t+1))/4.0E+0-(r*(1-s)*(t+1))/2.0E+0)/2.0E+0)-((1-s)*(t+1))/8.0E+0)+invJ.xelem(1,1)*((-((-((1-r)*(1-t_2))/4.0E+0)-((1-r)*s*(t+1))/2.0E+0-((1-r_2)*(t+1))/4.0E+0)/2.0E+0)-((1-r)*(t+1))/8.0E+0)+invJ.xelem(1,2)*(((1-r)*(1-s))/8.0E+0-((-((1-r)*(1-s)*t)/2.0E+0)+((1-r)*(1-s_2))/4.0E+0+((1-r_2)*(1-s))/4.0E+0)/2.0E+0);
          Bt.xelem(2,2) = invJ.xelem(2,0)*((-((-((1-s)*(1-t_2))/4.0E+0)-((1-s_2)*(t+1))/4.0E+0-(r*(1-s)*(t+1))/2.0E+0)/2.0E+0)-((1-s)*(t+1))/8.0E+0)+invJ.xelem(2,1)*((-((-((1-r)*(1-t_2))/4.0E+0)-((1-r)*s*(t+1))/2.0E+0-((1-r_2)*(t+1))/4.0E+0)/2.0E+0)-((1-r)*(t+1))/8.0E+0)+invJ.xelem(2,2)*(((1-r)*(1-s))/8.0E+0-((-((1-r)*(1-s)*t)/2.0E+0)+((1-r)*(1-s_2))/4.0E+0+((1-r_2)*(1-s))/4.0E+0)/2.0E+0);
          Bt.xelem(0,3) = invJ.xelem(0,0)*(((1-s)*(t+1))/8.0E+0-(((1-s)*(1-t_2))/4.0E+0+((1-s_2)*(t+1))/4.0E+0-(r*(1-s)*(t+1))/2.0E+0)/2.0E+0)+invJ.xelem(0,1)*((-((-((r+1)*(1-t_2))/4.0E+0)-((r+1)*s*(t+1))/2.0E+0-((1-r_2)*(t+1))/4.0E+0)/2.0E+0)-((r+1)*(t+1))/8.0E+0)+invJ.xelem(0,2)*(((r+1)*(1-s))/8.0E+0-((-((r+1)*(1-s)*t)/2.0E+0)+((r+1)*(1-s_2))/4.0E+0+((1-r_2)*(1-s))/4.0E+0)/2.0E+0);
          Bt.xelem(1,3) = invJ.xelem(1,0)*(((1-s)*(t+1))/8.0E+0-(((1-s)*(1-t_2))/4.0E+0+((1-s_2)*(t+1))/4.0E+0-(r*(1-s)*(t+1))/2.0E+0)/2.0E+0)+invJ.xelem(1,1)*((-((-((r+1)*(1-t_2))/4.0E+0)-((r+1)*s*(t+1))/2.0E+0-((1-r_2)*(t+1))/4.0E+0)/2.0E+0)-((r+1)*(t+1))/8.0E+0)+invJ.xelem(1,2)*(((r+1)*(1-s))/8.0E+0-((-((r+1)*(1-s)*t)/2.0E+0)+((r+1)*(1-s_2))/4.0E+0+((1-r_2)*(1-s))/4.0E+0)/2.0E+0);
          Bt.xelem(2,3) = invJ.xelem(2,0)*(((1-s)*(t+1))/8.0E+0-(((1-s)*(1-t_2))/4.0E+0+((1-s_2)*(t+1))/4.0E+0-(r*(1-s)*(t+1))/2.0E+0)/2.0E+0)+invJ.xelem(2,1)*((-((-((r+1)*(1-t_2))/4.0E+0)-((r+1)*s*(t+1))/2.0E+0-((1-r_2)*(t+1))/4.0E+0)/2.0E+0)-((r+1)*(t+1))/8.0E+0)+invJ.xelem(2,2)*(((r+1)*(1-s))/8.0E+0-((-((r+1)*(1-s)*t)/2.0E+0)+((r+1)*(1-s_2))/4.0E+0+((1-r_2)*(1-s))/4.0E+0)/2.0E+0);
          Bt.xelem(0,4) = invJ.xelem(0,0)*(((s+1)*(1-t))/8.0E+0-(((s+1)*(1-t_2))/4.0E+0+((1-s_2)*(1-t))/4.0E+0-(r*(s+1)*(1-t))/2.0E+0)/2.0E+0)+invJ.xelem(0,1)*(((r+1)*(1-t))/8.0E+0-(((r+1)*(1-t_2))/4.0E+0-((r+1)*s*(1-t))/2.0E+0+((1-r_2)*(1-t))/4.0E+0)/2.0E+0)+invJ.xelem(0,2)*((-((-((r+1)*(s+1)*t)/2.0E+0)-((r+1)*(1-s_2))/4.0E+0-((1-r_2)*(s+1))/4.0E+0)/2.0E+0)-((r+1)*(s+1))/8.0E+0);
          Bt.xelem(1,4) = invJ.xelem(1,0)*(((s+1)*(1-t))/8.0E+0-(((s+1)*(1-t_2))/4.0E+0+((1-s_2)*(1-t))/4.0E+0-(r*(s+1)*(1-t))/2.0E+0)/2.0E+0)+invJ.xelem(1,1)*(((r+1)*(1-t))/8.0E+0-(((r+1)*(1-t_2))/4.0E+0-((r+1)*s*(1-t))/2.0E+0+((1-r_2)*(1-t))/4.0E+0)/2.0E+0)+invJ.xelem(1,2)*((-((-((r+1)*(s+1)*t)/2.0E+0)-((r+1)*(1-s_2))/4.0E+0-((1-r_2)*(s+1))/4.0E+0)/2.0E+0)-((r+1)*(s+1))/8.0E+0);
          Bt.xelem(2,4) = invJ.xelem(2,0)*(((s+1)*(1-t))/8.0E+0-(((s+1)*(1-t_2))/4.0E+0+((1-s_2)*(1-t))/4.0E+0-(r*(s+1)*(1-t))/2.0E+0)/2.0E+0)+invJ.xelem(2,1)*(((r+1)*(1-t))/8.0E+0-(((r+1)*(1-t_2))/4.0E+0-((r+1)*s*(1-t))/2.0E+0+((1-r_2)*(1-t))/4.0E+0)/2.0E+0)+invJ.xelem(2,2)*((-((-((r+1)*(s+1)*t)/2.0E+0)-((r+1)*(1-s_2))/4.0E+0-((1-r_2)*(s+1))/4.0E+0)/2.0E+0)-((r+1)*(s+1))/8.0E+0);
          Bt.xelem(0,5) = invJ.xelem(0,0)*((-((-((s+1)*(1-t_2))/4.0E+0)-((1-s_2)*(1-t))/4.0E+0-(r*(s+1)*(1-t))/2.0E+0)/2.0E+0)-((s+1)*(1-t))/8.0E+0)+invJ.xelem(0,1)*(((1-r)*(1-t))/8.0E+0-(((1-r)*(1-t_2))/4.0E+0-((1-r)*s*(1-t))/2.0E+0+((1-r_2)*(1-t))/4.0E+0)/2.0E+0)+invJ.xelem(0,2)*((-((-((1-r)*(s+1)*t)/2.0E+0)-((1-r)*(1-s_2))/4.0E+0-((1-r_2)*(s+1))/4.0E+0)/2.0E+0)-((1-r)*(s+1))/8.0E+0);
          Bt.xelem(1,5) = invJ.xelem(1,0)*((-((-((s+1)*(1-t_2))/4.0E+0)-((1-s_2)*(1-t))/4.0E+0-(r*(s+1)*(1-t))/2.0E+0)/2.0E+0)-((s+1)*(1-t))/8.0E+0)+invJ.xelem(1,1)*(((1-r)*(1-t))/8.0E+0-(((1-r)*(1-t_2))/4.0E+0-((1-r)*s*(1-t))/2.0E+0+((1-r_2)*(1-t))/4.0E+0)/2.0E+0)+invJ.xelem(1,2)*((-((-((1-r)*(s+1)*t)/2.0E+0)-((1-r)*(1-s_2))/4.0E+0-((1-r_2)*(s+1))/4.0E+0)/2.0E+0)-((1-r)*(s+1))/8.0E+0);
          Bt.xelem(2,5) = invJ.xelem(2,0)*((-((-((s+1)*(1-t_2))/4.0E+0)-((1-s_2)*(1-t))/4.0E+0-(r*(s+1)*(1-t))/2.0E+0)/2.0E+0)-((s+1)*(1-t))/8.0E+0)+invJ.xelem(2,1)*(((1-r)*(1-t))/8.0E+0-(((1-r)*(1-t_2))/4.0E+0-((1-r)*s*(1-t))/2.0E+0+((1-r_2)*(1-t))/4.0E+0)/2.0E+0)+invJ.xelem(2,2)*((-((-((1-r)*(s+1)*t)/2.0E+0)-((1-r)*(1-s_2))/4.0E+0-((1-r_2)*(s+1))/4.0E+0)/2.0E+0)-((1-r)*(s+1))/8.0E+0);
          Bt.xelem(0,6) = invJ.xelem(0,0)*((-((-((1-s)*(1-t_2))/4.0E+0)-((1-s_2)*(1-t))/4.0E+0-(r*(1-s)*(1-t))/2.0E+0)/2.0E+0)-((1-s)*(1-t))/8.0E+0)+invJ.xelem(0,1)*((-((-((1-r)*(1-t_2))/4.0E+0)-((1-r)*s*(1-t))/2.0E+0-((1-r_2)*(1-t))/4.0E+0)/2.0E+0)-((1-r)*(1-t))/8.0E+0)+invJ.xelem(0,2)*((-((-((1-r)*(1-s)*t)/2.0E+0)-((1-r)*(1-s_2))/4.0E+0-((1-r_2)*(1-s))/4.0E+0)/2.0E+0)-((1-r)*(1-s))/8.0E+0);
          Bt.xelem(1,6) = invJ.xelem(1,0)*((-((-((1-s)*(1-t_2))/4.0E+0)-((1-s_2)*(1-t))/4.0E+0-(r*(1-s)*(1-t))/2.0E+0)/2.0E+0)-((1-s)*(1-t))/8.0E+0)+invJ.xelem(1,1)*((-((-((1-r)*(1-t_2))/4.0E+0)-((1-r)*s*(1-t))/2.0E+0-((1-r_2)*(1-t))/4.0E+0)/2.0E+0)-((1-r)*(1-t))/8.0E+0)+invJ.xelem(1,2)*((-((-((1-r)*(1-s)*t)/2.0E+0)-((1-r)*(1-s_2))/4.0E+0-((1-r_2)*(1-s))/4.0E+0)/2.0E+0)-((1-r)*(1-s))/8.0E+0);
          Bt.xelem(2,6) = invJ.xelem(2,0)*((-((-((1-s)*(1-t_2))/4.0E+0)-((1-s_2)*(1-t))/4.0E+0-(r*(1-s)*(1-t))/2.0E+0)/2.0E+0)-((1-s)*(1-t))/8.0E+0)+invJ.xelem(2,1)*((-((-((1-r)*(1-t_2))/4.0E+0)-((1-r)*s*(1-t))/2.0E+0-((1-r_2)*(1-t))/4.0E+0)/2.0E+0)-((1-r)*(1-t))/8.0E+0)+invJ.xelem(2,2)*((-((-((1-r)*(1-s)*t)/2.0E+0)-((1-r)*(1-s_2))/4.0E+0-((1-r_2)*(1-s))/4.0E+0)/2.0E+0)-((1-r)*(1-s))/8.0E+0);
          Bt.xelem(0,7) = invJ.xelem(0,0)*(((1-s)*(1-t))/8.0E+0-(((1-s)*(1-t_2))/4.0E+0+((1-s_2)*(1-t))/4.0E+0-(r*(1-s)*(1-t))/2.0E+0)/2.0E+0)+invJ.xelem(0,1)*((-((-((r+1)*(1-t_2))/4.0E+0)-((r+1)*s*(1-t))/2.0E+0-((1-r_2)*(1-t))/4.0E+0)/2.0E+0)-((r+1)*(1-t))/8.0E+0)+invJ.xelem(0,2)*((-((-((r+1)*(1-s)*t)/2.0E+0)-((r+1)*(1-s_2))/4.0E+0-((1-r_2)*(1-s))/4.0E+0)/2.0E+0)-((r+1)*(1-s))/8.0E+0);
          Bt.xelem(1,7) = invJ.xelem(1,0)*(((1-s)*(1-t))/8.0E+0-(((1-s)*(1-t_2))/4.0E+0+((1-s_2)*(1-t))/4.0E+0-(r*(1-s)*(1-t))/2.0E+0)/2.0E+0)+invJ.xelem(1,1)*((-((-((r+1)*(1-t_2))/4.0E+0)-((r+1)*s*(1-t))/2.0E+0-((1-r_2)*(1-t))/4.0E+0)/2.0E+0)-((r+1)*(1-t))/8.0E+0)+invJ.xelem(1,2)*((-((-((r+1)*(1-s)*t)/2.0E+0)-((r+1)*(1-s_2))/4.0E+0-((1-r_2)*(1-s))/4.0E+0)/2.0E+0)-((r+1)*(1-s))/8.0E+0);
          Bt.xelem(2,7) = invJ.xelem(2,0)*(((1-s)*(1-t))/8.0E+0-(((1-s)*(1-t_2))/4.0E+0+((1-s_2)*(1-t))/4.0E+0-(r*(1-s)*(1-t))/2.0E+0)/2.0E+0)+invJ.xelem(2,1)*((-((-((r+1)*(1-t_2))/4.0E+0)-((r+1)*s*(1-t))/2.0E+0-((1-r_2)*(1-t))/4.0E+0)/2.0E+0)-((r+1)*(1-t))/8.0E+0)+invJ.xelem(2,2)*((-((-((r+1)*(1-s)*t)/2.0E+0)-((r+1)*(1-s_2))/4.0E+0-((1-r_2)*(1-s))/4.0E+0)/2.0E+0)-((r+1)*(1-s))/8.0E+0);
          Bt.xelem(0,8) = (-(invJ.xelem(0,0)*r*(s+1)*(t+1))/2.0E+0)+(invJ.xelem(0,1)*(1-r_2)*(t+1))/4.0E+0+(invJ.xelem(0,2)*(1-r_2)*(s+1))/4.0E+0;
          Bt.xelem(1,8) = (-(invJ.xelem(1,0)*r*(s+1)*(t+1))/2.0E+0)+(invJ.xelem(1,1)*(1-r_2)*(t+1))/4.0E+0+(invJ.xelem(1,2)*(1-r_2)*(s+1))/4.0E+0;
          Bt.xelem(2,8) = (-(invJ.xelem(2,0)*r*(s+1)*(t+1))/2.0E+0)+(invJ.xelem(2,1)*(1-r_2)*(t+1))/4.0E+0+(invJ.xelem(2,2)*(1-r_2)*(s+1))/4.0E+0;
          Bt.xelem(0,9) = (-(invJ.xelem(0,0)*(1-s_2)*(t+1))/4.0E+0)-(invJ.xelem(0,1)*(1-r)*s*(t+1))/2.0E+0+(invJ.xelem(0,2)*(1-r)*(1-s_2))/4.0E+0;
          Bt.xelem(1,9) = (-(invJ.xelem(1,0)*(1-s_2)*(t+1))/4.0E+0)-(invJ.xelem(1,1)*(1-r)*s*(t+1))/2.0E+0+(invJ.xelem(1,2)*(1-r)*(1-s_2))/4.0E+0;
          Bt.xelem(2,9) = (-(invJ.xelem(2,0)*(1-s_2)*(t+1))/4.0E+0)-(invJ.xelem(2,1)*(1-r)*s*(t+1))/2.0E+0+(invJ.xelem(2,2)*(1-r)*(1-s_2))/4.0E+0;
          Bt.xelem(0,10) = (-(invJ.xelem(0,0)*r*(1-s)*(t+1))/2.0E+0)-(invJ.xelem(0,1)*(1-r_2)*(t+1))/4.0E+0+(invJ.xelem(0,2)*(1-r_2)*(1-s))/4.0E+0;
          Bt.xelem(1,10) = (-(invJ.xelem(1,0)*r*(1-s)*(t+1))/2.0E+0)-(invJ.xelem(1,1)*(1-r_2)*(t+1))/4.0E+0+(invJ.xelem(1,2)*(1-r_2)*(1-s))/4.0E+0;
          Bt.xelem(2,10) = (-(invJ.xelem(2,0)*r*(1-s)*(t+1))/2.0E+0)-(invJ.xelem(2,1)*(1-r_2)*(t+1))/4.0E+0+(invJ.xelem(2,2)*(1-r_2)*(1-s))/4.0E+0;
          Bt.xelem(0,11) = (invJ.xelem(0,0)*(1-s_2)*(t+1))/4.0E+0-(invJ.xelem(0,1)*(r+1)*s*(t+1))/2.0E+0+(invJ.xelem(0,2)*(r+1)*(1-s_2))/4.0E+0;
          Bt.xelem(1,11) = (invJ.xelem(1,0)*(1-s_2)*(t+1))/4.0E+0-(invJ.xelem(1,1)*(r+1)*s*(t+1))/2.0E+0+(invJ.xelem(1,2)*(r+1)*(1-s_2))/4.0E+0;
          Bt.xelem(2,11) = (invJ.xelem(2,0)*(1-s_2)*(t+1))/4.0E+0-(invJ.xelem(2,1)*(r+1)*s*(t+1))/2.0E+0+(invJ.xelem(2,2)*(r+1)*(1-s_2))/4.0E+0;
          Bt.xelem(0,12) = (-(invJ.xelem(0,0)*r*(s+1)*(1-t))/2.0E+0)+(invJ.xelem(0,1)*(1-r_2)*(1-t))/4.0E+0-(invJ.xelem(0,2)*(1-r_2)*(s+1))/4.0E+0;
          Bt.xelem(1,12) = (-(invJ.xelem(1,0)*r*(s+1)*(1-t))/2.0E+0)+(invJ.xelem(1,1)*(1-r_2)*(1-t))/4.0E+0-(invJ.xelem(1,2)*(1-r_2)*(s+1))/4.0E+0;
          Bt.xelem(2,12) = (-(invJ.xelem(2,0)*r*(s+1)*(1-t))/2.0E+0)+(invJ.xelem(2,1)*(1-r_2)*(1-t))/4.0E+0-(invJ.xelem(2,2)*(1-r_2)*(s+1))/4.0E+0;
          Bt.xelem(0,13) = (-(invJ.xelem(0,0)*(1-s_2)*(1-t))/4.0E+0)-(invJ.xelem(0,1)*(1-r)*s*(1-t))/2.0E+0-(invJ.xelem(0,2)*(1-r)*(1-s_2))/4.0E+0;
          Bt.xelem(1,13) = (-(invJ.xelem(1,0)*(1-s_2)*(1-t))/4.0E+0)-(invJ.xelem(1,1)*(1-r)*s*(1-t))/2.0E+0-(invJ.xelem(1,2)*(1-r)*(1-s_2))/4.0E+0;
          Bt.xelem(2,13) = (-(invJ.xelem(2,0)*(1-s_2)*(1-t))/4.0E+0)-(invJ.xelem(2,1)*(1-r)*s*(1-t))/2.0E+0-(invJ.xelem(2,2)*(1-r)*(1-s_2))/4.0E+0;
          Bt.xelem(0,14) = (-(invJ.xelem(0,0)*r*(1-s)*(1-t))/2.0E+0)-(invJ.xelem(0,1)*(1-r_2)*(1-t))/4.0E+0-(invJ.xelem(0,2)*(1-r_2)*(1-s))/4.0E+0;
          Bt.xelem(1,14) = (-(invJ.xelem(1,0)*r*(1-s)*(1-t))/2.0E+0)-(invJ.xelem(1,1)*(1-r_2)*(1-t))/4.0E+0-(invJ.xelem(1,2)*(1-r_2)*(1-s))/4.0E+0;
          Bt.xelem(2,14) = (-(invJ.xelem(2,0)*r*(1-s)*(1-t))/2.0E+0)-(invJ.xelem(2,1)*(1-r_2)*(1-t))/4.0E+0-(invJ.xelem(2,2)*(1-r_2)*(1-s))/4.0E+0;
          Bt.xelem(0,15) = (invJ.xelem(0,0)*(1-s_2)*(1-t))/4.0E+0-(invJ.xelem(0,1)*(r+1)*s*(1-t))/2.0E+0-(invJ.xelem(0,2)*(r+1)*(1-s_2))/4.0E+0;
          Bt.xelem(1,15) = (invJ.xelem(1,0)*(1-s_2)*(1-t))/4.0E+0-(invJ.xelem(1,1)*(r+1)*s*(1-t))/2.0E+0-(invJ.xelem(1,2)*(r+1)*(1-s_2))/4.0E+0;
          Bt.xelem(2,15) = (invJ.xelem(2,0)*(1-s_2)*(1-t))/4.0E+0-(invJ.xelem(2,1)*(r+1)*s*(1-t))/2.0E+0-(invJ.xelem(2,2)*(r+1)*(1-s_2))/4.0E+0;
          Bt.xelem(0,16) = (invJ.xelem(0,0)*(s+1)*(1-t_2))/4.0E+0+(invJ.xelem(0,1)*(r+1)*(1-t_2))/4.0E+0-(invJ.xelem(0,2)*(r+1)*(s+1)*t)/2.0E+0;
          Bt.xelem(1,16) = (invJ.xelem(1,0)*(s+1)*(1-t_2))/4.0E+0+(invJ.xelem(1,1)*(r+1)*(1-t_2))/4.0E+0-(invJ.xelem(1,2)*(r+1)*(s+1)*t)/2.0E+0;
          Bt.xelem(2,16) = (invJ.xelem(2,0)*(s+1)*(1-t_2))/4.0E+0+(invJ.xelem(2,1)*(r+1)*(1-t_2))/4.0E+0-(invJ.xelem(2,2)*(r+1)*(s+1)*t)/2.0E+0;
          Bt.xelem(0,17) = (-(invJ.xelem(0,0)*(s+1)*(1-t_2))/4.0E+0)+(invJ.xelem(0,1)*(1-r)*(1-t_2))/4.0E+0-(invJ.xelem(0,2)*(1-r)*(s+1)*t)/2.0E+0;
          Bt.xelem(1,17) = (-(invJ.xelem(1,0)*(s+1)*(1-t_2))/4.0E+0)+(invJ.xelem(1,1)*(1-r)*(1-t_2))/4.0E+0-(invJ.xelem(1,2)*(1-r)*(s+1)*t)/2.0E+0;
          Bt.xelem(2,17) = (-(invJ.xelem(2,0)*(s+1)*(1-t_2))/4.0E+0)+(invJ.xelem(2,1)*(1-r)*(1-t_2))/4.0E+0-(invJ.xelem(2,2)*(1-r)*(s+1)*t)/2.0E+0;
          Bt.xelem(0,18) = (-(invJ.xelem(0,0)*(1-s)*(1-t_2))/4.0E+0)-(invJ.xelem(0,1)*(1-r)*(1-t_2))/4.0E+0-(invJ.xelem(0,2)*(1-r)*(1-s)*t)/2.0E+0;
          Bt.xelem(1,18) = (-(invJ.xelem(1,0)*(1-s)*(1-t_2))/4.0E+0)-(invJ.xelem(1,1)*(1-r)*(1-t_2))/4.0E+0-(invJ.xelem(1,2)*(1-r)*(1-s)*t)/2.0E+0;
          Bt.xelem(2,18) = (-(invJ.xelem(2,0)*(1-s)*(1-t_2))/4.0E+0)-(invJ.xelem(2,1)*(1-r)*(1-t_2))/4.0E+0-(invJ.xelem(2,2)*(1-r)*(1-s)*t)/2.0E+0;
          Bt.xelem(0,19) = (invJ.xelem(0,0)*(1-s)*(1-t_2))/4.0E+0-(invJ.xelem(0,1)*(r+1)*(1-t_2))/4.0E+0-(invJ.xelem(0,2)*(r+1)*(1-s)*t)/2.0E+0;
          Bt.xelem(1,19) = (invJ.xelem(1,0)*(1-s)*(1-t_2))/4.0E+0-(invJ.xelem(1,1)*(r+1)*(1-t_2))/4.0E+0-(invJ.xelem(1,2)*(r+1)*(1-s)*t)/2.0E+0;
          Bt.xelem(2,19) = (invJ.xelem(2,0)*(1-s)*(1-t_2))/4.0E+0-(invJ.xelem(2,1)*(r+1)*(1-t_2))/4.0E+0-(invJ.xelem(2,2)*(r+1)*(1-s)*t)/2.0E+0;
     }     
};

class Penta15: public Element3D
{
public:
     Penta15(ElementTypes::TypeId eltype, octave_idx_type id, const Matrix& X, const Material* material, const int32NDArray& nodes, const StrainField& oRefStrain)
          :Element3D(eltype, id, X, material, nodes, oRefStrain) {
          FEM_ASSERT(nodes.numel() == 15);
     }

     virtual const IntegrationRule& GetIntegrationRule(FemMatrixType eMatType) const {
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
          static constexpr double w[2][M] = {{w1, w2, w3, w4, w5, w6, w7}, {1./3., 1./3., 1./3., 1./3., 1./3., 1./3.}};
          static const double t[2][N] = {{0.774596669241483, 0., -0.774596669241483}, {1., 0., -1.}};
          static const double alpha[2][N] = {{0.555555555555556, 0.888888888888889, 0.555555555555556}, {2./3., 2./3., 2./3.}};

          static array<IntegrationRule, 2> rgIntegRule;

          octave_idx_type iIntegRule, iNumPoints;

          switch (eMatType) {
          case MAT_MASS_LUMPED:
               iIntegRule = 1;
               iNumPoints = 6;
               break;
          default:
               iIntegRule = 0;
               iNumPoints = M;
          }

          if (!rgIntegRule[iIntegRule].iGetNumEvalPoints()) {
               rgIntegRule[iIntegRule].SetNumEvalPoints(M * N, 3);

               for (octave_idx_type i = 0; i < iNumPoints; ++i) {
                    for (octave_idx_type j = 0; j < N; ++j) {
                         rgIntegRule[iIntegRule].SetWeight(i * N + j, 0.5 * w[iIntegRule][i] * alpha[iIntegRule][j]);
                         rgIntegRule[iIntegRule].SetPosition(i * N + j, 0, r[iIntegRule][i]);
                         rgIntegRule[iIntegRule].SetPosition(i * N + j, 1, s[iIntegRule][i]);
                         rgIntegRule[iIntegRule].SetPosition(i * N + j, 2, t[iIntegRule][j]);
                    }
               }
          }

          return rgIntegRule[iIntegRule];
     }

protected:
     virtual double Jacobian(const ColumnVector& rv, Matrix& J) const {
          FEM_ASSERT(J.rows() == 3);
          FEM_ASSERT(J.columns() == 3);
          FEM_ASSERT(rv.numel() == 3);

          const double r = rv.xelem(0);
          const double s = rv.xelem(1);
          const double t = rv.xelem(2);
          const double t2 = t * t;

          J.xelem(0,0) = X.xelem(0,13)*(1-t2)-X.xelem(0,12)*(1-t2)-(X.xelem(0,3)*(t+1)*(t-2*s-2*r))/2.0E+0+(X.xelem(0,4)*(t+1)*(t+2*r-2))/2.0E+0-2*X.xelem(0,11)*s*(t+1)+2*X.xelem(0,10)*s*(t+1)+2*X.xelem(0,9)*((-s)-r+1)*(t+1)-X.xelem(0,3)*((-s)-r+1)*(t+1)-2*X.xelem(0,9)*r*(t+1)+X.xelem(0,4)*r*(t+1)-(X.xelem(0,0)*(1-t)*((-t)-2*s-2*r))/2.0E+0+(X.xelem(0,1)*(1-t)*((-t)+2*r-2))/2.0E+0-2*X.xelem(0,8)*s*(1-t)+2*X.xelem(0,7)*s*(1-t)+2*X.xelem(0,6)*((-s)-r+1)*(1-t)-X.xelem(0,0)*((-s)-r+1)*(1-t)-2*X.xelem(0,6)*r*(1-t)+X.xelem(0,1)*r*(1-t);
          J.xelem(1,0) = X.xelem(0,14)*(1-t2)-X.xelem(0,12)*(1-t2)+(X.xelem(0,5)*(t+1)*(t+2*s-2))/2.0E+0-(X.xelem(0,3)*(t+1)*(t-2*s-2*r))/2.0E+0-2*X.xelem(0,11)*s*(t+1)+X.xelem(0,5)*s*(t+1)+2*X.xelem(0,11)*((-s)-r+1)*(t+1)-X.xelem(0,3)*((-s)-r+1)*(t+1)+2*X.xelem(0,10)*r*(t+1)-2*X.xelem(0,9)*r*(t+1)+(X.xelem(0,2)*(1-t)*((-t)+2*s-2))/2.0E+0-(X.xelem(0,0)*(1-t)*((-t)-2*s-2*r))/2.0E+0-2*X.xelem(0,8)*s*(1-t)+X.xelem(0,2)*s*(1-t)+2*X.xelem(0,8)*((-s)-r+1)*(1-t)-X.xelem(0,0)*((-s)-r+1)*(1-t)+2*X.xelem(0,7)*r*(1-t)-2*X.xelem(0,6)*r*(1-t);
          J.xelem(2,0) = (X.xelem(0,5)*s*(t+2*s-2))/2.0E+0+(X.xelem(0,3)*((-s)-r+1)*(t-2*s-2*r))/2.0E+0+(X.xelem(0,4)*r*(t+2*r-2))/2.0E+0+(X.xelem(0,5)*s*(t+1))/2.0E+0+(X.xelem(0,3)*((-s)-r+1)*(t+1))/2.0E+0+(X.xelem(0,4)*r*(t+1))/2.0E+0-2*X.xelem(0,14)*s*t-2*X.xelem(0,12)*((-s)-r+1)*t-2*X.xelem(0,13)*r*t-(X.xelem(0,2)*s*((-t)+2*s-2))/2.0E+0-(X.xelem(0,0)*((-s)-r+1)*((-t)-2*s-2*r))/2.0E+0-(X.xelem(0,1)*r*((-t)+2*r-2))/2.0E+0-(X.xelem(0,2)*s*(1-t))/2.0E+0-(X.xelem(0,0)*((-s)-r+1)*(1-t))/2.0E+0-(X.xelem(0,1)*r*(1-t))/2.0E+0+2*X.xelem(0,11)*((-s)-r+1)*s-2*X.xelem(0,8)*((-s)-r+1)*s+2*X.xelem(0,10)*r*s-2*X.xelem(0,7)*r*s+2*X.xelem(0,9)*r*((-s)-r+1)-2*X.xelem(0,6)*r*((-s)-r+1);
          J.xelem(0,1) = X.xelem(1,13)*(1-t2)-X.xelem(1,12)*(1-t2)-(X.xelem(1,3)*(t+1)*(t-2*s-2*r))/2.0E+0+(X.xelem(1,4)*(t+1)*(t+2*r-2))/2.0E+0-2*X.xelem(1,11)*s*(t+1)+2*X.xelem(1,10)*s*(t+1)+2*X.xelem(1,9)*((-s)-r+1)*(t+1)-X.xelem(1,3)*((-s)-r+1)*(t+1)-2*X.xelem(1,9)*r*(t+1)+X.xelem(1,4)*r*(t+1)-(X.xelem(1,0)*(1-t)*((-t)-2*s-2*r))/2.0E+0+(X.xelem(1,1)*(1-t)*((-t)+2*r-2))/2.0E+0-2*X.xelem(1,8)*s*(1-t)+2*X.xelem(1,7)*s*(1-t)+2*X.xelem(1,6)*((-s)-r+1)*(1-t)-X.xelem(1,0)*((-s)-r+1)*(1-t)-2*X.xelem(1,6)*r*(1-t)+X.xelem(1,1)*r*(1-t);
          J.xelem(1,1) = X.xelem(1,14)*(1-t2)-X.xelem(1,12)*(1-t2)+(X.xelem(1,5)*(t+1)*(t+2*s-2))/2.0E+0-(X.xelem(1,3)*(t+1)*(t-2*s-2*r))/2.0E+0-2*X.xelem(1,11)*s*(t+1)+X.xelem(1,5)*s*(t+1)+2*X.xelem(1,11)*((-s)-r+1)*(t+1)-X.xelem(1,3)*((-s)-r+1)*(t+1)+2*X.xelem(1,10)*r*(t+1)-2*X.xelem(1,9)*r*(t+1)+(X.xelem(1,2)*(1-t)*((-t)+2*s-2))/2.0E+0-(X.xelem(1,0)*(1-t)*((-t)-2*s-2*r))/2.0E+0-2*X.xelem(1,8)*s*(1-t)+X.xelem(1,2)*s*(1-t)+2*X.xelem(1,8)*((-s)-r+1)*(1-t)-X.xelem(1,0)*((-s)-r+1)*(1-t)+2*X.xelem(1,7)*r*(1-t)-2*X.xelem(1,6)*r*(1-t);
          J.xelem(2,1) = (X.xelem(1,5)*s*(t+2*s-2))/2.0E+0+(X.xelem(1,3)*((-s)-r+1)*(t-2*s-2*r))/2.0E+0+(X.xelem(1,4)*r*(t+2*r-2))/2.0E+0+(X.xelem(1,5)*s*(t+1))/2.0E+0+(X.xelem(1,3)*((-s)-r+1)*(t+1))/2.0E+0+(X.xelem(1,4)*r*(t+1))/2.0E+0-2*X.xelem(1,14)*s*t-2*X.xelem(1,12)*((-s)-r+1)*t-2*X.xelem(1,13)*r*t-(X.xelem(1,2)*s*((-t)+2*s-2))/2.0E+0-(X.xelem(1,0)*((-s)-r+1)*((-t)-2*s-2*r))/2.0E+0-(X.xelem(1,1)*r*((-t)+2*r-2))/2.0E+0-(X.xelem(1,2)*s*(1-t))/2.0E+0-(X.xelem(1,0)*((-s)-r+1)*(1-t))/2.0E+0-(X.xelem(1,1)*r*(1-t))/2.0E+0+2*X.xelem(1,11)*((-s)-r+1)*s-2*X.xelem(1,8)*((-s)-r+1)*s+2*X.xelem(1,10)*r*s-2*X.xelem(1,7)*r*s+2*X.xelem(1,9)*r*((-s)-r+1)-2*X.xelem(1,6)*r*((-s)-r+1);
          J.xelem(0,2) = X.xelem(2,13)*(1-t2)-X.xelem(2,12)*(1-t2)-(X.xelem(2,3)*(t+1)*(t-2*s-2*r))/2.0E+0+(X.xelem(2,4)*(t+1)*(t+2*r-2))/2.0E+0-2*X.xelem(2,11)*s*(t+1)+2*X.xelem(2,10)*s*(t+1)+2*X.xelem(2,9)*((-s)-r+1)*(t+1)-X.xelem(2,3)*((-s)-r+1)*(t+1)-2*X.xelem(2,9)*r*(t+1)+X.xelem(2,4)*r*(t+1)-(X.xelem(2,0)*(1-t)*((-t)-2*s-2*r))/2.0E+0+(X.xelem(2,1)*(1-t)*((-t)+2*r-2))/2.0E+0-2*X.xelem(2,8)*s*(1-t)+2*X.xelem(2,7)*s*(1-t)+2*X.xelem(2,6)*((-s)-r+1)*(1-t)-X.xelem(2,0)*((-s)-r+1)*(1-t)-2*X.xelem(2,6)*r*(1-t)+X.xelem(2,1)*r*(1-t);
          J.xelem(1,2) = X.xelem(2,14)*(1-t2)-X.xelem(2,12)*(1-t2)+(X.xelem(2,5)*(t+1)*(t+2*s-2))/2.0E+0-(X.xelem(2,3)*(t+1)*(t-2*s-2*r))/2.0E+0-2*X.xelem(2,11)*s*(t+1)+X.xelem(2,5)*s*(t+1)+2*X.xelem(2,11)*((-s)-r+1)*(t+1)-X.xelem(2,3)*((-s)-r+1)*(t+1)+2*X.xelem(2,10)*r*(t+1)-2*X.xelem(2,9)*r*(t+1)+(X.xelem(2,2)*(1-t)*((-t)+2*s-2))/2.0E+0-(X.xelem(2,0)*(1-t)*((-t)-2*s-2*r))/2.0E+0-2*X.xelem(2,8)*s*(1-t)+X.xelem(2,2)*s*(1-t)+2*X.xelem(2,8)*((-s)-r+1)*(1-t)-X.xelem(2,0)*((-s)-r+1)*(1-t)+2*X.xelem(2,7)*r*(1-t)-2*X.xelem(2,6)*r*(1-t);
          J.xelem(2,2) = (X.xelem(2,5)*s*(t+2*s-2))/2.0E+0+(X.xelem(2,3)*((-s)-r+1)*(t-2*s-2*r))/2.0E+0+(X.xelem(2,4)*r*(t+2*r-2))/2.0E+0+(X.xelem(2,5)*s*(t+1))/2.0E+0+(X.xelem(2,3)*((-s)-r+1)*(t+1))/2.0E+0+(X.xelem(2,4)*r*(t+1))/2.0E+0-2*X.xelem(2,14)*s*t-2*X.xelem(2,12)*((-s)-r+1)*t-2*X.xelem(2,13)*r*t-(X.xelem(2,2)*s*((-t)+2*s-2))/2.0E+0-(X.xelem(2,0)*((-s)-r+1)*((-t)-2*s-2*r))/2.0E+0-(X.xelem(2,1)*r*((-t)+2*r-2))/2.0E+0-(X.xelem(2,2)*s*(1-t))/2.0E+0-(X.xelem(2,0)*((-s)-r+1)*(1-t))/2.0E+0-(X.xelem(2,1)*r*(1-t))/2.0E+0+2*X.xelem(2,11)*((-s)-r+1)*s-2*X.xelem(2,8)*((-s)-r+1)*s+2*X.xelem(2,10)*r*s-2*X.xelem(2,7)*r*s+2*X.xelem(2,9)*r*((-s)-r+1)-2*X.xelem(2,6)*r*((-s)-r+1);

          return J.determinant();
     }

     virtual void DispInterpMatrix(const ColumnVector& rv, Matrix& H) const {
          FEM_ASSERT(rv.numel() == 3);
          FEM_ASSERT(H.rows() == 3);
          FEM_ASSERT(H.columns() == 45);

          const double r = rv.xelem(0);
          const double s = rv.xelem(1);
          const double t = rv.xelem(2);
          const double t2 = t * t;

          H.xelem(0,0) = (((-s)-r+1)*(1-t)*((-t)-2*s-2*r))/2.0E+0;
          H.xelem(1,0) = 0;
          H.xelem(2,0) = 0;
          H.xelem(0,1) = 0;
          H.xelem(1,1) = (((-s)-r+1)*(1-t)*((-t)-2*s-2*r))/2.0E+0;
          H.xelem(2,1) = 0;
          H.xelem(0,2) = 0;
          H.xelem(1,2) = 0;
          H.xelem(2,2) = (((-s)-r+1)*(1-t)*((-t)-2*s-2*r))/2.0E+0;
          H.xelem(0,3) = (r*(1-t)*((-t)+2*r-2))/2.0E+0;
          H.xelem(1,3) = 0;
          H.xelem(2,3) = 0;
          H.xelem(0,4) = 0;
          H.xelem(1,4) = (r*(1-t)*((-t)+2*r-2))/2.0E+0;
          H.xelem(2,4) = 0;
          H.xelem(0,5) = 0;
          H.xelem(1,5) = 0;
          H.xelem(2,5) = (r*(1-t)*((-t)+2*r-2))/2.0E+0;
          H.xelem(0,6) = (s*(1-t)*((-t)+2*s-2))/2.0E+0;
          H.xelem(1,6) = 0;
          H.xelem(2,6) = 0;
          H.xelem(0,7) = 0;
          H.xelem(1,7) = (s*(1-t)*((-t)+2*s-2))/2.0E+0;
          H.xelem(2,7) = 0;
          H.xelem(0,8) = 0;
          H.xelem(1,8) = 0;
          H.xelem(2,8) = (s*(1-t)*((-t)+2*s-2))/2.0E+0;
          H.xelem(0,9) = (((-s)-r+1)*(t+1)*(t-2*s-2*r))/2.0E+0;
          H.xelem(1,9) = 0;
          H.xelem(2,9) = 0;
          H.xelem(0,10) = 0;
          H.xelem(1,10) = (((-s)-r+1)*(t+1)*(t-2*s-2*r))/2.0E+0;
          H.xelem(2,10) = 0;
          H.xelem(0,11) = 0;
          H.xelem(1,11) = 0;
          H.xelem(2,11) = (((-s)-r+1)*(t+1)*(t-2*s-2*r))/2.0E+0;
          H.xelem(0,12) = (r*(t+1)*(t+2*r-2))/2.0E+0;
          H.xelem(1,12) = 0;
          H.xelem(2,12) = 0;
          H.xelem(0,13) = 0;
          H.xelem(1,13) = (r*(t+1)*(t+2*r-2))/2.0E+0;
          H.xelem(2,13) = 0;
          H.xelem(0,14) = 0;
          H.xelem(1,14) = 0;
          H.xelem(2,14) = (r*(t+1)*(t+2*r-2))/2.0E+0;
          H.xelem(0,15) = (s*(t+1)*(t+2*s-2))/2.0E+0;
          H.xelem(1,15) = 0;
          H.xelem(2,15) = 0;
          H.xelem(0,16) = 0;
          H.xelem(1,16) = (s*(t+1)*(t+2*s-2))/2.0E+0;
          H.xelem(2,16) = 0;
          H.xelem(0,17) = 0;
          H.xelem(1,17) = 0;
          H.xelem(2,17) = (s*(t+1)*(t+2*s-2))/2.0E+0;
          H.xelem(0,18) = 2*r*((-s)-r+1)*(1-t);
          H.xelem(1,18) = 0;
          H.xelem(2,18) = 0;
          H.xelem(0,19) = 0;
          H.xelem(1,19) = 2*r*((-s)-r+1)*(1-t);
          H.xelem(2,19) = 0;
          H.xelem(0,20) = 0;
          H.xelem(1,20) = 0;
          H.xelem(2,20) = 2*r*((-s)-r+1)*(1-t);
          H.xelem(0,21) = 2*r*s*(1-t);
          H.xelem(1,21) = 0;
          H.xelem(2,21) = 0;
          H.xelem(0,22) = 0;
          H.xelem(1,22) = 2*r*s*(1-t);
          H.xelem(2,22) = 0;
          H.xelem(0,23) = 0;
          H.xelem(1,23) = 0;
          H.xelem(2,23) = 2*r*s*(1-t);
          H.xelem(0,24) = 2*((-s)-r+1)*s*(1-t);
          H.xelem(1,24) = 0;
          H.xelem(2,24) = 0;
          H.xelem(0,25) = 0;
          H.xelem(1,25) = 2*((-s)-r+1)*s*(1-t);
          H.xelem(2,25) = 0;
          H.xelem(0,26) = 0;
          H.xelem(1,26) = 0;
          H.xelem(2,26) = 2*((-s)-r+1)*s*(1-t);
          H.xelem(0,27) = 2*r*((-s)-r+1)*(t+1);
          H.xelem(1,27) = 0;
          H.xelem(2,27) = 0;
          H.xelem(0,28) = 0;
          H.xelem(1,28) = 2*r*((-s)-r+1)*(t+1);
          H.xelem(2,28) = 0;
          H.xelem(0,29) = 0;
          H.xelem(1,29) = 0;
          H.xelem(2,29) = 2*r*((-s)-r+1)*(t+1);
          H.xelem(0,30) = 2*r*s*(t+1);
          H.xelem(1,30) = 0;
          H.xelem(2,30) = 0;
          H.xelem(0,31) = 0;
          H.xelem(1,31) = 2*r*s*(t+1);
          H.xelem(2,31) = 0;
          H.xelem(0,32) = 0;
          H.xelem(1,32) = 0;
          H.xelem(2,32) = 2*r*s*(t+1);
          H.xelem(0,33) = 2*((-s)-r+1)*s*(t+1);
          H.xelem(1,33) = 0;
          H.xelem(2,33) = 0;
          H.xelem(0,34) = 0;
          H.xelem(1,34) = 2*((-s)-r+1)*s*(t+1);
          H.xelem(2,34) = 0;
          H.xelem(0,35) = 0;
          H.xelem(1,35) = 0;
          H.xelem(2,35) = 2*((-s)-r+1)*s*(t+1);
          H.xelem(0,36) = ((-s)-r+1)*(1-t2);
          H.xelem(1,36) = 0;
          H.xelem(2,36) = 0;
          H.xelem(0,37) = 0;
          H.xelem(1,37) = ((-s)-r+1)*(1-t2);
          H.xelem(2,37) = 0;
          H.xelem(0,38) = 0;
          H.xelem(1,38) = 0;
          H.xelem(2,38) = ((-s)-r+1)*(1-t2);
          H.xelem(0,39) = r*(1-t2);
          H.xelem(1,39) = 0;
          H.xelem(2,39) = 0;
          H.xelem(0,40) = 0;
          H.xelem(1,40) = r*(1-t2);
          H.xelem(2,40) = 0;
          H.xelem(0,41) = 0;
          H.xelem(1,41) = 0;
          H.xelem(2,41) = r*(1-t2);
          H.xelem(0,42) = s*(1-t2);
          H.xelem(1,42) = 0;
          H.xelem(2,42) = 0;
          H.xelem(0,43) = 0;
          H.xelem(1,43) = s*(1-t2);
          H.xelem(2,43) = 0;
          H.xelem(0,44) = 0;
          H.xelem(1,44) = 0;
          H.xelem(2,44) = s*(1-t2);
     }

     virtual void StrainMatrix(const ColumnVector& rv, const Matrix& invJ, Matrix& B) const {
          FEM_ASSERT(rv.numel() == 3);
          FEM_ASSERT(invJ.rows() == 3);
          FEM_ASSERT(invJ.columns() == 3);
          FEM_ASSERT(B.rows() == 6);
          FEM_ASSERT(B.columns() == 45);

          const double r = rv.xelem(0);
          const double s = rv.xelem(1);
          const double t = rv.xelem(2);
          const double t2 = t * t;

          B.xelem(0,0) = invJ(0,1)*((-((1-t)*((-t)-2*s-2*r))/2.0E+0)-((-s)-r+1)*(1-t))+invJ(0,0)*((-((1-t)*((-t)-2*s-2*r))/2.0E+0)-((-s)-r+1)*(1-t))+invJ(0,2)*((-(((-s)-r+1)*((-t)-2*s-2*r))/2.0E+0)-(((-s)-r+1)*(1-t))/2.0E+0);
          B.xelem(1,0) = 0;
          B.xelem(2,0) = 0;
          B.xelem(3,0) = invJ(1,1)*((-((1-t)*((-t)-2*s-2*r))/2.0E+0)-((-s)-r+1)*(1-t))+invJ(1,0)*((-((1-t)*((-t)-2*s-2*r))/2.0E+0)-((-s)-r+1)*(1-t))+invJ(1,2)*((-(((-s)-r+1)*((-t)-2*s-2*r))/2.0E+0)-(((-s)-r+1)*(1-t))/2.0E+0);
          B.xelem(4,0) = 0;
          B.xelem(5,0) = invJ(2,1)*((-((1-t)*((-t)-2*s-2*r))/2.0E+0)-((-s)-r+1)*(1-t))+invJ(2,0)*((-((1-t)*((-t)-2*s-2*r))/2.0E+0)-((-s)-r+1)*(1-t))+invJ(2,2)*((-(((-s)-r+1)*((-t)-2*s-2*r))/2.0E+0)-(((-s)-r+1)*(1-t))/2.0E+0);
          B.xelem(0,1) = 0;
          B.xelem(1,1) = invJ(1,1)*((-((1-t)*((-t)-2*s-2*r))/2.0E+0)-((-s)-r+1)*(1-t))+invJ(1,0)*((-((1-t)*((-t)-2*s-2*r))/2.0E+0)-((-s)-r+1)*(1-t))+invJ(1,2)*((-(((-s)-r+1)*((-t)-2*s-2*r))/2.0E+0)-(((-s)-r+1)*(1-t))/2.0E+0);
          B.xelem(2,1) = 0;
          B.xelem(3,1) = invJ(0,1)*((-((1-t)*((-t)-2*s-2*r))/2.0E+0)-((-s)-r+1)*(1-t))+invJ(0,0)*((-((1-t)*((-t)-2*s-2*r))/2.0E+0)-((-s)-r+1)*(1-t))+invJ(0,2)*((-(((-s)-r+1)*((-t)-2*s-2*r))/2.0E+0)-(((-s)-r+1)*(1-t))/2.0E+0);
          B.xelem(4,1) = invJ(2,1)*((-((1-t)*((-t)-2*s-2*r))/2.0E+0)-((-s)-r+1)*(1-t))+invJ(2,0)*((-((1-t)*((-t)-2*s-2*r))/2.0E+0)-((-s)-r+1)*(1-t))+invJ(2,2)*((-(((-s)-r+1)*((-t)-2*s-2*r))/2.0E+0)-(((-s)-r+1)*(1-t))/2.0E+0);
          B.xelem(5,1) = 0;
          B.xelem(0,2) = 0;
          B.xelem(1,2) = 0;
          B.xelem(2,2) = invJ(2,1)*((-((1-t)*((-t)-2*s-2*r))/2.0E+0)-((-s)-r+1)*(1-t))+invJ(2,0)*((-((1-t)*((-t)-2*s-2*r))/2.0E+0)-((-s)-r+1)*(1-t))+invJ(2,2)*((-(((-s)-r+1)*((-t)-2*s-2*r))/2.0E+0)-(((-s)-r+1)*(1-t))/2.0E+0);
          B.xelem(3,2) = 0;
          B.xelem(4,2) = invJ(1,1)*((-((1-t)*((-t)-2*s-2*r))/2.0E+0)-((-s)-r+1)*(1-t))+invJ(1,0)*((-((1-t)*((-t)-2*s-2*r))/2.0E+0)-((-s)-r+1)*(1-t))+invJ(1,2)*((-(((-s)-r+1)*((-t)-2*s-2*r))/2.0E+0)-(((-s)-r+1)*(1-t))/2.0E+0);
          B.xelem(5,2) = invJ(0,1)*((-((1-t)*((-t)-2*s-2*r))/2.0E+0)-((-s)-r+1)*(1-t))+invJ(0,0)*((-((1-t)*((-t)-2*s-2*r))/2.0E+0)-((-s)-r+1)*(1-t))+invJ(0,2)*((-(((-s)-r+1)*((-t)-2*s-2*r))/2.0E+0)-(((-s)-r+1)*(1-t))/2.0E+0);
          B.xelem(0,3) = invJ(0,0)*(((1-t)*((-t)+2*r-2))/2.0E+0+r*(1-t))+invJ(0,2)*((-(r*((-t)+2*r-2))/2.0E+0)-(r*(1-t))/2.0E+0);
          B.xelem(1,3) = 0;
          B.xelem(2,3) = 0;
          B.xelem(3,3) = invJ(1,0)*(((1-t)*((-t)+2*r-2))/2.0E+0+r*(1-t))+invJ(1,2)*((-(r*((-t)+2*r-2))/2.0E+0)-(r*(1-t))/2.0E+0);
          B.xelem(4,3) = 0;
          B.xelem(5,3) = invJ(2,0)*(((1-t)*((-t)+2*r-2))/2.0E+0+r*(1-t))+invJ(2,2)*((-(r*((-t)+2*r-2))/2.0E+0)-(r*(1-t))/2.0E+0);
          B.xelem(0,4) = 0;
          B.xelem(1,4) = invJ(1,0)*(((1-t)*((-t)+2*r-2))/2.0E+0+r*(1-t))+invJ(1,2)*((-(r*((-t)+2*r-2))/2.0E+0)-(r*(1-t))/2.0E+0);
          B.xelem(2,4) = 0;
          B.xelem(3,4) = invJ(0,0)*(((1-t)*((-t)+2*r-2))/2.0E+0+r*(1-t))+invJ(0,2)*((-(r*((-t)+2*r-2))/2.0E+0)-(r*(1-t))/2.0E+0);
          B.xelem(4,4) = invJ(2,0)*(((1-t)*((-t)+2*r-2))/2.0E+0+r*(1-t))+invJ(2,2)*((-(r*((-t)+2*r-2))/2.0E+0)-(r*(1-t))/2.0E+0);
          B.xelem(5,4) = 0;
          B.xelem(0,5) = 0;
          B.xelem(1,5) = 0;
          B.xelem(2,5) = invJ(2,0)*(((1-t)*((-t)+2*r-2))/2.0E+0+r*(1-t))+invJ(2,2)*((-(r*((-t)+2*r-2))/2.0E+0)-(r*(1-t))/2.0E+0);
          B.xelem(3,5) = 0;
          B.xelem(4,5) = invJ(1,0)*(((1-t)*((-t)+2*r-2))/2.0E+0+r*(1-t))+invJ(1,2)*((-(r*((-t)+2*r-2))/2.0E+0)-(r*(1-t))/2.0E+0);
          B.xelem(5,5) = invJ(0,0)*(((1-t)*((-t)+2*r-2))/2.0E+0+r*(1-t))+invJ(0,2)*((-(r*((-t)+2*r-2))/2.0E+0)-(r*(1-t))/2.0E+0);
          B.xelem(0,6) = invJ(0,1)*(((1-t)*((-t)+2*s-2))/2.0E+0+s*(1-t))+invJ(0,2)*((-(s*((-t)+2*s-2))/2.0E+0)-(s*(1-t))/2.0E+0);
          B.xelem(1,6) = 0;
          B.xelem(2,6) = 0;
          B.xelem(3,6) = invJ(1,1)*(((1-t)*((-t)+2*s-2))/2.0E+0+s*(1-t))+invJ(1,2)*((-(s*((-t)+2*s-2))/2.0E+0)-(s*(1-t))/2.0E+0);
          B.xelem(4,6) = 0;
          B.xelem(5,6) = invJ(2,1)*(((1-t)*((-t)+2*s-2))/2.0E+0+s*(1-t))+invJ(2,2)*((-(s*((-t)+2*s-2))/2.0E+0)-(s*(1-t))/2.0E+0);
          B.xelem(0,7) = 0;
          B.xelem(1,7) = invJ(1,1)*(((1-t)*((-t)+2*s-2))/2.0E+0+s*(1-t))+invJ(1,2)*((-(s*((-t)+2*s-2))/2.0E+0)-(s*(1-t))/2.0E+0);
          B.xelem(2,7) = 0;
          B.xelem(3,7) = invJ(0,1)*(((1-t)*((-t)+2*s-2))/2.0E+0+s*(1-t))+invJ(0,2)*((-(s*((-t)+2*s-2))/2.0E+0)-(s*(1-t))/2.0E+0);
          B.xelem(4,7) = invJ(2,1)*(((1-t)*((-t)+2*s-2))/2.0E+0+s*(1-t))+invJ(2,2)*((-(s*((-t)+2*s-2))/2.0E+0)-(s*(1-t))/2.0E+0);
          B.xelem(5,7) = 0;
          B.xelem(0,8) = 0;
          B.xelem(1,8) = 0;
          B.xelem(2,8) = invJ(2,1)*(((1-t)*((-t)+2*s-2))/2.0E+0+s*(1-t))+invJ(2,2)*((-(s*((-t)+2*s-2))/2.0E+0)-(s*(1-t))/2.0E+0);
          B.xelem(3,8) = 0;
          B.xelem(4,8) = invJ(1,1)*(((1-t)*((-t)+2*s-2))/2.0E+0+s*(1-t))+invJ(1,2)*((-(s*((-t)+2*s-2))/2.0E+0)-(s*(1-t))/2.0E+0);
          B.xelem(5,8) = invJ(0,1)*(((1-t)*((-t)+2*s-2))/2.0E+0+s*(1-t))+invJ(0,2)*((-(s*((-t)+2*s-2))/2.0E+0)-(s*(1-t))/2.0E+0);
          B.xelem(0,9) = invJ(0,1)*((-((t+1)*(t-2*s-2*r))/2.0E+0)-((-s)-r+1)*(t+1))+invJ(0,0)*((-((t+1)*(t-2*s-2*r))/2.0E+0)-((-s)-r+1)*(t+1))+invJ(0,2)*((((-s)-r+1)*(t-2*s-2*r))/2.0E+0+(((-s)-r+1)*(t+1))/2.0E+0);
          B.xelem(1,9) = 0;
          B.xelem(2,9) = 0;
          B.xelem(3,9) = invJ(1,1)*((-((t+1)*(t-2*s-2*r))/2.0E+0)-((-s)-r+1)*(t+1))+invJ(1,0)*((-((t+1)*(t-2*s-2*r))/2.0E+0)-((-s)-r+1)*(t+1))+invJ(1,2)*((((-s)-r+1)*(t-2*s-2*r))/2.0E+0+(((-s)-r+1)*(t+1))/2.0E+0);
          B.xelem(4,9) = 0;
          B.xelem(5,9) = invJ(2,1)*((-((t+1)*(t-2*s-2*r))/2.0E+0)-((-s)-r+1)*(t+1))+invJ(2,0)*((-((t+1)*(t-2*s-2*r))/2.0E+0)-((-s)-r+1)*(t+1))+invJ(2,2)*((((-s)-r+1)*(t-2*s-2*r))/2.0E+0+(((-s)-r+1)*(t+1))/2.0E+0);
          B.xelem(0,10) = 0;
          B.xelem(1,10) = invJ(1,1)*((-((t+1)*(t-2*s-2*r))/2.0E+0)-((-s)-r+1)*(t+1))+invJ(1,0)*((-((t+1)*(t-2*s-2*r))/2.0E+0)-((-s)-r+1)*(t+1))+invJ(1,2)*((((-s)-r+1)*(t-2*s-2*r))/2.0E+0+(((-s)-r+1)*(t+1))/2.0E+0);
          B.xelem(2,10) = 0;
          B.xelem(3,10) = invJ(0,1)*((-((t+1)*(t-2*s-2*r))/2.0E+0)-((-s)-r+1)*(t+1))+invJ(0,0)*((-((t+1)*(t-2*s-2*r))/2.0E+0)-((-s)-r+1)*(t+1))+invJ(0,2)*((((-s)-r+1)*(t-2*s-2*r))/2.0E+0+(((-s)-r+1)*(t+1))/2.0E+0);
          B.xelem(4,10) = invJ(2,1)*((-((t+1)*(t-2*s-2*r))/2.0E+0)-((-s)-r+1)*(t+1))+invJ(2,0)*((-((t+1)*(t-2*s-2*r))/2.0E+0)-((-s)-r+1)*(t+1))+invJ(2,2)*((((-s)-r+1)*(t-2*s-2*r))/2.0E+0+(((-s)-r+1)*(t+1))/2.0E+0);
          B.xelem(5,10) = 0;
          B.xelem(0,11) = 0;
          B.xelem(1,11) = 0;
          B.xelem(2,11) = invJ(2,1)*((-((t+1)*(t-2*s-2*r))/2.0E+0)-((-s)-r+1)*(t+1))+invJ(2,0)*((-((t+1)*(t-2*s-2*r))/2.0E+0)-((-s)-r+1)*(t+1))+invJ(2,2)*((((-s)-r+1)*(t-2*s-2*r))/2.0E+0+(((-s)-r+1)*(t+1))/2.0E+0);
          B.xelem(3,11) = 0;
          B.xelem(4,11) = invJ(1,1)*((-((t+1)*(t-2*s-2*r))/2.0E+0)-((-s)-r+1)*(t+1))+invJ(1,0)*((-((t+1)*(t-2*s-2*r))/2.0E+0)-((-s)-r+1)*(t+1))+invJ(1,2)*((((-s)-r+1)*(t-2*s-2*r))/2.0E+0+(((-s)-r+1)*(t+1))/2.0E+0);
          B.xelem(5,11) = invJ(0,1)*((-((t+1)*(t-2*s-2*r))/2.0E+0)-((-s)-r+1)*(t+1))+invJ(0,0)*((-((t+1)*(t-2*s-2*r))/2.0E+0)-((-s)-r+1)*(t+1))+invJ(0,2)*((((-s)-r+1)*(t-2*s-2*r))/2.0E+0+(((-s)-r+1)*(t+1))/2.0E+0);
          B.xelem(0,12) = invJ(0,0)*(((t+1)*(t+2*r-2))/2.0E+0+r*(t+1))+invJ(0,2)*((r*(t+2*r-2))/2.0E+0+(r*(t+1))/2.0E+0);
          B.xelem(1,12) = 0;
          B.xelem(2,12) = 0;
          B.xelem(3,12) = invJ(1,0)*(((t+1)*(t+2*r-2))/2.0E+0+r*(t+1))+invJ(1,2)*((r*(t+2*r-2))/2.0E+0+(r*(t+1))/2.0E+0);
          B.xelem(4,12) = 0;
          B.xelem(5,12) = invJ(2,0)*(((t+1)*(t+2*r-2))/2.0E+0+r*(t+1))+invJ(2,2)*((r*(t+2*r-2))/2.0E+0+(r*(t+1))/2.0E+0);
          B.xelem(0,13) = 0;
          B.xelem(1,13) = invJ(1,0)*(((t+1)*(t+2*r-2))/2.0E+0+r*(t+1))+invJ(1,2)*((r*(t+2*r-2))/2.0E+0+(r*(t+1))/2.0E+0);
          B.xelem(2,13) = 0;
          B.xelem(3,13) = invJ(0,0)*(((t+1)*(t+2*r-2))/2.0E+0+r*(t+1))+invJ(0,2)*((r*(t+2*r-2))/2.0E+0+(r*(t+1))/2.0E+0);
          B.xelem(4,13) = invJ(2,0)*(((t+1)*(t+2*r-2))/2.0E+0+r*(t+1))+invJ(2,2)*((r*(t+2*r-2))/2.0E+0+(r*(t+1))/2.0E+0);
          B.xelem(5,13) = 0;
          B.xelem(0,14) = 0;
          B.xelem(1,14) = 0;
          B.xelem(2,14) = invJ(2,0)*(((t+1)*(t+2*r-2))/2.0E+0+r*(t+1))+invJ(2,2)*((r*(t+2*r-2))/2.0E+0+(r*(t+1))/2.0E+0);
          B.xelem(3,14) = 0;
          B.xelem(4,14) = invJ(1,0)*(((t+1)*(t+2*r-2))/2.0E+0+r*(t+1))+invJ(1,2)*((r*(t+2*r-2))/2.0E+0+(r*(t+1))/2.0E+0);
          B.xelem(5,14) = invJ(0,0)*(((t+1)*(t+2*r-2))/2.0E+0+r*(t+1))+invJ(0,2)*((r*(t+2*r-2))/2.0E+0+(r*(t+1))/2.0E+0);
          B.xelem(0,15) = invJ(0,1)*(((t+1)*(t+2*s-2))/2.0E+0+s*(t+1))+invJ(0,2)*((s*(t+2*s-2))/2.0E+0+(s*(t+1))/2.0E+0);
          B.xelem(1,15) = 0;
          B.xelem(2,15) = 0;
          B.xelem(3,15) = invJ(1,1)*(((t+1)*(t+2*s-2))/2.0E+0+s*(t+1))+invJ(1,2)*((s*(t+2*s-2))/2.0E+0+(s*(t+1))/2.0E+0);
          B.xelem(4,15) = 0;
          B.xelem(5,15) = invJ(2,1)*(((t+1)*(t+2*s-2))/2.0E+0+s*(t+1))+invJ(2,2)*((s*(t+2*s-2))/2.0E+0+(s*(t+1))/2.0E+0);
          B.xelem(0,16) = 0;
          B.xelem(1,16) = invJ(1,1)*(((t+1)*(t+2*s-2))/2.0E+0+s*(t+1))+invJ(1,2)*((s*(t+2*s-2))/2.0E+0+(s*(t+1))/2.0E+0);
          B.xelem(2,16) = 0;
          B.xelem(3,16) = invJ(0,1)*(((t+1)*(t+2*s-2))/2.0E+0+s*(t+1))+invJ(0,2)*((s*(t+2*s-2))/2.0E+0+(s*(t+1))/2.0E+0);
          B.xelem(4,16) = invJ(2,1)*(((t+1)*(t+2*s-2))/2.0E+0+s*(t+1))+invJ(2,2)*((s*(t+2*s-2))/2.0E+0+(s*(t+1))/2.0E+0);
          B.xelem(5,16) = 0;
          B.xelem(0,17) = 0;
          B.xelem(1,17) = 0;
          B.xelem(2,17) = invJ(2,1)*(((t+1)*(t+2*s-2))/2.0E+0+s*(t+1))+invJ(2,2)*((s*(t+2*s-2))/2.0E+0+(s*(t+1))/2.0E+0);
          B.xelem(3,17) = 0;
          B.xelem(4,17) = invJ(1,1)*(((t+1)*(t+2*s-2))/2.0E+0+s*(t+1))+invJ(1,2)*((s*(t+2*s-2))/2.0E+0+(s*(t+1))/2.0E+0);
          B.xelem(5,17) = invJ(0,1)*(((t+1)*(t+2*s-2))/2.0E+0+s*(t+1))+invJ(0,2)*((s*(t+2*s-2))/2.0E+0+(s*(t+1))/2.0E+0);
          B.xelem(0,18) = (-2*invJ(0,1)*r*(1-t))+invJ(0,0)*(2*((-s)-r+1)*(1-t)-2*r*(1-t))-2*invJ(0,2)*r*((-s)-r+1);
          B.xelem(1,18) = 0;
          B.xelem(2,18) = 0;
          B.xelem(3,18) = (-2*invJ(1,1)*r*(1-t))+invJ(1,0)*(2*((-s)-r+1)*(1-t)-2*r*(1-t))-2*invJ(1,2)*r*((-s)-r+1);
          B.xelem(4,18) = 0;
          B.xelem(5,18) = (-2*invJ(2,1)*r*(1-t))+invJ(2,0)*(2*((-s)-r+1)*(1-t)-2*r*(1-t))-2*invJ(2,2)*r*((-s)-r+1);
          B.xelem(0,19) = 0;
          B.xelem(1,19) = (-2*invJ(1,1)*r*(1-t))+invJ(1,0)*(2*((-s)-r+1)*(1-t)-2*r*(1-t))-2*invJ(1,2)*r*((-s)-r+1);
          B.xelem(2,19) = 0;
          B.xelem(3,19) = (-2*invJ(0,1)*r*(1-t))+invJ(0,0)*(2*((-s)-r+1)*(1-t)-2*r*(1-t))-2*invJ(0,2)*r*((-s)-r+1);
          B.xelem(4,19) = (-2*invJ(2,1)*r*(1-t))+invJ(2,0)*(2*((-s)-r+1)*(1-t)-2*r*(1-t))-2*invJ(2,2)*r*((-s)-r+1);
          B.xelem(5,19) = 0;
          B.xelem(0,20) = 0;
          B.xelem(1,20) = 0;
          B.xelem(2,20) = (-2*invJ(2,1)*r*(1-t))+invJ(2,0)*(2*((-s)-r+1)*(1-t)-2*r*(1-t))-2*invJ(2,2)*r*((-s)-r+1);
          B.xelem(3,20) = 0;
          B.xelem(4,20) = (-2*invJ(1,1)*r*(1-t))+invJ(1,0)*(2*((-s)-r+1)*(1-t)-2*r*(1-t))-2*invJ(1,2)*r*((-s)-r+1);
          B.xelem(5,20) = (-2*invJ(0,1)*r*(1-t))+invJ(0,0)*(2*((-s)-r+1)*(1-t)-2*r*(1-t))-2*invJ(0,2)*r*((-s)-r+1);
          B.xelem(0,21) = 2*invJ(0,0)*s*(1-t)+2*invJ(0,1)*r*(1-t)-2*invJ(0,2)*r*s;
          B.xelem(1,21) = 0;
          B.xelem(2,21) = 0;
          B.xelem(3,21) = 2*invJ(1,0)*s*(1-t)+2*invJ(1,1)*r*(1-t)-2*invJ(1,2)*r*s;
          B.xelem(4,21) = 0;
          B.xelem(5,21) = 2*invJ(2,0)*s*(1-t)+2*invJ(2,1)*r*(1-t)-2*invJ(2,2)*r*s;
          B.xelem(0,22) = 0;
          B.xelem(1,22) = 2*invJ(1,0)*s*(1-t)+2*invJ(1,1)*r*(1-t)-2*invJ(1,2)*r*s;
          B.xelem(2,22) = 0;
          B.xelem(3,22) = 2*invJ(0,0)*s*(1-t)+2*invJ(0,1)*r*(1-t)-2*invJ(0,2)*r*s;
          B.xelem(4,22) = 2*invJ(2,0)*s*(1-t)+2*invJ(2,1)*r*(1-t)-2*invJ(2,2)*r*s;
          B.xelem(5,22) = 0;
          B.xelem(0,23) = 0;
          B.xelem(1,23) = 0;
          B.xelem(2,23) = 2*invJ(2,0)*s*(1-t)+2*invJ(2,1)*r*(1-t)-2*invJ(2,2)*r*s;
          B.xelem(3,23) = 0;
          B.xelem(4,23) = 2*invJ(1,0)*s*(1-t)+2*invJ(1,1)*r*(1-t)-2*invJ(1,2)*r*s;
          B.xelem(5,23) = 2*invJ(0,0)*s*(1-t)+2*invJ(0,1)*r*(1-t)-2*invJ(0,2)*r*s;
          B.xelem(0,24) = (-2*invJ(0,0)*s*(1-t))+invJ(0,1)*(2*((-s)-r+1)*(1-t)-2*s*(1-t))-2*invJ(0,2)*((-s)-r+1)*s;
          B.xelem(1,24) = 0;
          B.xelem(2,24) = 0;
          B.xelem(3,24) = (-2*invJ(1,0)*s*(1-t))+invJ(1,1)*(2*((-s)-r+1)*(1-t)-2*s*(1-t))-2*invJ(1,2)*((-s)-r+1)*s;
          B.xelem(4,24) = 0;
          B.xelem(5,24) = (-2*invJ(2,0)*s*(1-t))+invJ(2,1)*(2*((-s)-r+1)*(1-t)-2*s*(1-t))-2*invJ(2,2)*((-s)-r+1)*s;
          B.xelem(0,25) = 0;
          B.xelem(1,25) = (-2*invJ(1,0)*s*(1-t))+invJ(1,1)*(2*((-s)-r+1)*(1-t)-2*s*(1-t))-2*invJ(1,2)*((-s)-r+1)*s;
          B.xelem(2,25) = 0;
          B.xelem(3,25) = (-2*invJ(0,0)*s*(1-t))+invJ(0,1)*(2*((-s)-r+1)*(1-t)-2*s*(1-t))-2*invJ(0,2)*((-s)-r+1)*s;
          B.xelem(4,25) = (-2*invJ(2,0)*s*(1-t))+invJ(2,1)*(2*((-s)-r+1)*(1-t)-2*s*(1-t))-2*invJ(2,2)*((-s)-r+1)*s;
          B.xelem(5,25) = 0;
          B.xelem(0,26) = 0;
          B.xelem(1,26) = 0;
          B.xelem(2,26) = (-2*invJ(2,0)*s*(1-t))+invJ(2,1)*(2*((-s)-r+1)*(1-t)-2*s*(1-t))-2*invJ(2,2)*((-s)-r+1)*s;
          B.xelem(3,26) = 0;
          B.xelem(4,26) = (-2*invJ(1,0)*s*(1-t))+invJ(1,1)*(2*((-s)-r+1)*(1-t)-2*s*(1-t))-2*invJ(1,2)*((-s)-r+1)*s;
          B.xelem(5,26) = (-2*invJ(0,0)*s*(1-t))+invJ(0,1)*(2*((-s)-r+1)*(1-t)-2*s*(1-t))-2*invJ(0,2)*((-s)-r+1)*s;
          B.xelem(0,27) = invJ(0,0)*(2*((-s)-r+1)*(t+1)-2*r*(t+1))-2*invJ(0,1)*r*(t+1)+2*invJ(0,2)*r*((-s)-r+1);
          B.xelem(1,27) = 0;
          B.xelem(2,27) = 0;
          B.xelem(3,27) = invJ(1,0)*(2*((-s)-r+1)*(t+1)-2*r*(t+1))-2*invJ(1,1)*r*(t+1)+2*invJ(1,2)*r*((-s)-r+1);
          B.xelem(4,27) = 0;
          B.xelem(5,27) = invJ(2,0)*(2*((-s)-r+1)*(t+1)-2*r*(t+1))-2*invJ(2,1)*r*(t+1)+2*invJ(2,2)*r*((-s)-r+1);
          B.xelem(0,28) = 0;
          B.xelem(1,28) = invJ(1,0)*(2*((-s)-r+1)*(t+1)-2*r*(t+1))-2*invJ(1,1)*r*(t+1)+2*invJ(1,2)*r*((-s)-r+1);
          B.xelem(2,28) = 0;
          B.xelem(3,28) = invJ(0,0)*(2*((-s)-r+1)*(t+1)-2*r*(t+1))-2*invJ(0,1)*r*(t+1)+2*invJ(0,2)*r*((-s)-r+1);
          B.xelem(4,28) = invJ(2,0)*(2*((-s)-r+1)*(t+1)-2*r*(t+1))-2*invJ(2,1)*r*(t+1)+2*invJ(2,2)*r*((-s)-r+1);
          B.xelem(5,28) = 0;
          B.xelem(0,29) = 0;
          B.xelem(1,29) = 0;
          B.xelem(2,29) = invJ(2,0)*(2*((-s)-r+1)*(t+1)-2*r*(t+1))-2*invJ(2,1)*r*(t+1)+2*invJ(2,2)*r*((-s)-r+1);
          B.xelem(3,29) = 0;
          B.xelem(4,29) = invJ(1,0)*(2*((-s)-r+1)*(t+1)-2*r*(t+1))-2*invJ(1,1)*r*(t+1)+2*invJ(1,2)*r*((-s)-r+1);
          B.xelem(5,29) = invJ(0,0)*(2*((-s)-r+1)*(t+1)-2*r*(t+1))-2*invJ(0,1)*r*(t+1)+2*invJ(0,2)*r*((-s)-r+1);
          B.xelem(0,30) = 2*invJ(0,0)*s*(t+1)+2*invJ(0,1)*r*(t+1)+2*invJ(0,2)*r*s;
          B.xelem(1,30) = 0;
          B.xelem(2,30) = 0;
          B.xelem(3,30) = 2*invJ(1,0)*s*(t+1)+2*invJ(1,1)*r*(t+1)+2*invJ(1,2)*r*s;
          B.xelem(4,30) = 0;
          B.xelem(5,30) = 2*invJ(2,0)*s*(t+1)+2*invJ(2,1)*r*(t+1)+2*invJ(2,2)*r*s;
          B.xelem(0,31) = 0;
          B.xelem(1,31) = 2*invJ(1,0)*s*(t+1)+2*invJ(1,1)*r*(t+1)+2*invJ(1,2)*r*s;
          B.xelem(2,31) = 0;
          B.xelem(3,31) = 2*invJ(0,0)*s*(t+1)+2*invJ(0,1)*r*(t+1)+2*invJ(0,2)*r*s;
          B.xelem(4,31) = 2*invJ(2,0)*s*(t+1)+2*invJ(2,1)*r*(t+1)+2*invJ(2,2)*r*s;
          B.xelem(5,31) = 0;
          B.xelem(0,32) = 0;
          B.xelem(1,32) = 0;
          B.xelem(2,32) = 2*invJ(2,0)*s*(t+1)+2*invJ(2,1)*r*(t+1)+2*invJ(2,2)*r*s;
          B.xelem(3,32) = 0;
          B.xelem(4,32) = 2*invJ(1,0)*s*(t+1)+2*invJ(1,1)*r*(t+1)+2*invJ(1,2)*r*s;
          B.xelem(5,32) = 2*invJ(0,0)*s*(t+1)+2*invJ(0,1)*r*(t+1)+2*invJ(0,2)*r*s;
          B.xelem(0,33) = invJ(0,1)*(2*((-s)-r+1)*(t+1)-2*s*(t+1))-2*invJ(0,0)*s*(t+1)+2*invJ(0,2)*((-s)-r+1)*s;
          B.xelem(1,33) = 0;
          B.xelem(2,33) = 0;
          B.xelem(3,33) = invJ(1,1)*(2*((-s)-r+1)*(t+1)-2*s*(t+1))-2*invJ(1,0)*s*(t+1)+2*invJ(1,2)*((-s)-r+1)*s;
          B.xelem(4,33) = 0;
          B.xelem(5,33) = invJ(2,1)*(2*((-s)-r+1)*(t+1)-2*s*(t+1))-2*invJ(2,0)*s*(t+1)+2*invJ(2,2)*((-s)-r+1)*s;
          B.xelem(0,34) = 0;
          B.xelem(1,34) = invJ(1,1)*(2*((-s)-r+1)*(t+1)-2*s*(t+1))-2*invJ(1,0)*s*(t+1)+2*invJ(1,2)*((-s)-r+1)*s;
          B.xelem(2,34) = 0;
          B.xelem(3,34) = invJ(0,1)*(2*((-s)-r+1)*(t+1)-2*s*(t+1))-2*invJ(0,0)*s*(t+1)+2*invJ(0,2)*((-s)-r+1)*s;
          B.xelem(4,34) = invJ(2,1)*(2*((-s)-r+1)*(t+1)-2*s*(t+1))-2*invJ(2,0)*s*(t+1)+2*invJ(2,2)*((-s)-r+1)*s;
          B.xelem(5,34) = 0;
          B.xelem(0,35) = 0;
          B.xelem(1,35) = 0;
          B.xelem(2,35) = invJ(2,1)*(2*((-s)-r+1)*(t+1)-2*s*(t+1))-2*invJ(2,0)*s*(t+1)+2*invJ(2,2)*((-s)-r+1)*s;
          B.xelem(3,35) = 0;
          B.xelem(4,35) = invJ(1,1)*(2*((-s)-r+1)*(t+1)-2*s*(t+1))-2*invJ(1,0)*s*(t+1)+2*invJ(1,2)*((-s)-r+1)*s;
          B.xelem(5,35) = invJ(0,1)*(2*((-s)-r+1)*(t+1)-2*s*(t+1))-2*invJ(0,0)*s*(t+1)+2*invJ(0,2)*((-s)-r+1)*s;
          B.xelem(0,36) = invJ(0,1)*(t2-1)+invJ(0,0)*(t2-1)-2*invJ(0,2)*((-s)-r+1)*t;
          B.xelem(1,36) = 0;
          B.xelem(2,36) = 0;
          B.xelem(3,36) = invJ(1,1)*(t2-1)+invJ(1,0)*(t2-1)-2*invJ(1,2)*((-s)-r+1)*t;
          B.xelem(4,36) = 0;
          B.xelem(5,36) = invJ(2,1)*(t2-1)+invJ(2,0)*(t2-1)-2*invJ(2,2)*((-s)-r+1)*t;
          B.xelem(0,37) = 0;
          B.xelem(1,37) = invJ(1,1)*(t2-1)+invJ(1,0)*(t2-1)-2*invJ(1,2)*((-s)-r+1)*t;
          B.xelem(2,37) = 0;
          B.xelem(3,37) = invJ(0,1)*(t2-1)+invJ(0,0)*(t2-1)-2*invJ(0,2)*((-s)-r+1)*t;
          B.xelem(4,37) = invJ(2,1)*(t2-1)+invJ(2,0)*(t2-1)-2*invJ(2,2)*((-s)-r+1)*t;
          B.xelem(5,37) = 0;
          B.xelem(0,38) = 0;
          B.xelem(1,38) = 0;
          B.xelem(2,38) = invJ(2,1)*(t2-1)+invJ(2,0)*(t2-1)-2*invJ(2,2)*((-s)-r+1)*t;
          B.xelem(3,38) = 0;
          B.xelem(4,38) = invJ(1,1)*(t2-1)+invJ(1,0)*(t2-1)-2*invJ(1,2)*((-s)-r+1)*t;
          B.xelem(5,38) = invJ(0,1)*(t2-1)+invJ(0,0)*(t2-1)-2*invJ(0,2)*((-s)-r+1)*t;
          B.xelem(0,39) = invJ(0,0)*(1-t2)-2*invJ(0,2)*r*t;
          B.xelem(1,39) = 0;
          B.xelem(2,39) = 0;
          B.xelem(3,39) = invJ(1,0)*(1-t2)-2*invJ(1,2)*r*t;
          B.xelem(4,39) = 0;
          B.xelem(5,39) = invJ(2,0)*(1-t2)-2*invJ(2,2)*r*t;
          B.xelem(0,40) = 0;
          B.xelem(1,40) = invJ(1,0)*(1-t2)-2*invJ(1,2)*r*t;
          B.xelem(2,40) = 0;
          B.xelem(3,40) = invJ(0,0)*(1-t2)-2*invJ(0,2)*r*t;
          B.xelem(4,40) = invJ(2,0)*(1-t2)-2*invJ(2,2)*r*t;
          B.xelem(5,40) = 0;
          B.xelem(0,41) = 0;
          B.xelem(1,41) = 0;
          B.xelem(2,41) = invJ(2,0)*(1-t2)-2*invJ(2,2)*r*t;
          B.xelem(3,41) = 0;
          B.xelem(4,41) = invJ(1,0)*(1-t2)-2*invJ(1,2)*r*t;
          B.xelem(5,41) = invJ(0,0)*(1-t2)-2*invJ(0,2)*r*t;
          B.xelem(0,42) = invJ(0,1)*(1-t2)-2*invJ(0,2)*s*t;
          B.xelem(1,42) = 0;
          B.xelem(2,42) = 0;
          B.xelem(3,42) = invJ(1,1)*(1-t2)-2*invJ(1,2)*s*t;
          B.xelem(4,42) = 0;
          B.xelem(5,42) = invJ(2,1)*(1-t2)-2*invJ(2,2)*s*t;
          B.xelem(0,43) = 0;
          B.xelem(1,43) = invJ(1,1)*(1-t2)-2*invJ(1,2)*s*t;
          B.xelem(2,43) = 0;
          B.xelem(3,43) = invJ(0,1)*(1-t2)-2*invJ(0,2)*s*t;
          B.xelem(4,43) = invJ(2,1)*(1-t2)-2*invJ(2,2)*s*t;
          B.xelem(5,43) = 0;
          B.xelem(0,44) = 0;
          B.xelem(1,44) = 0;
          B.xelem(2,44) = invJ(2,1)*(1-t2)-2*invJ(2,2)*s*t;
          B.xelem(3,44) = 0;
          B.xelem(4,44) = invJ(1,1)*(1-t2)-2*invJ(1,2)*s*t;
          B.xelem(5,44) = invJ(0,1)*(1-t2)-2*invJ(0,2)*s*t;
     }

     virtual void TemperatureGradientMatrix(const ColumnVector& rv, const Matrix& invJ, Matrix& Bt) const {
          FEM_ASSERT(Bt.rows() == 3);
          FEM_ASSERT(Bt.columns() == 15);
          FEM_ASSERT(invJ.rows() == 3);
          FEM_ASSERT(invJ.columns() == 3);

          const double r = rv.xelem(0);
          const double s = rv.xelem(1);
          const double t = rv.xelem(2);
          const double t_2 = t * t;
          
          Bt.xelem(0,0) = invJ.xelem(0,1)*((-((1-t)*((-t)-2*s-2*r))/2.0E+0)-((-s)-r+1)*(1-t))+invJ.xelem(0,0)*((-((1-t)*((-t)-2*s-2*r))/2.0E+0)-((-s)-r+1)*(1-t))+invJ.xelem(0,2)*((-(((-s)-r+1)*((-t)-2*s-2*r))/2.0E+0)-(((-s)-r+1)*(1-t))/2.0E+0);
          Bt.xelem(1,0) = invJ.xelem(1,1)*((-((1-t)*((-t)-2*s-2*r))/2.0E+0)-((-s)-r+1)*(1-t))+invJ.xelem(1,0)*((-((1-t)*((-t)-2*s-2*r))/2.0E+0)-((-s)-r+1)*(1-t))+invJ.xelem(1,2)*((-(((-s)-r+1)*((-t)-2*s-2*r))/2.0E+0)-(((-s)-r+1)*(1-t))/2.0E+0);
          Bt.xelem(2,0) = invJ.xelem(2,1)*((-((1-t)*((-t)-2*s-2*r))/2.0E+0)-((-s)-r+1)*(1-t))+invJ.xelem(2,0)*((-((1-t)*((-t)-2*s-2*r))/2.0E+0)-((-s)-r+1)*(1-t))+invJ.xelem(2,2)*((-(((-s)-r+1)*((-t)-2*s-2*r))/2.0E+0)-(((-s)-r+1)*(1-t))/2.0E+0);
          Bt.xelem(0,1) = invJ.xelem(0,0)*(((1-t)*((-t)+2*r-2))/2.0E+0+r*(1-t))+invJ.xelem(0,2)*((-(r*((-t)+2*r-2))/2.0E+0)-(r*(1-t))/2.0E+0);
          Bt.xelem(1,1) = invJ.xelem(1,0)*(((1-t)*((-t)+2*r-2))/2.0E+0+r*(1-t))+invJ.xelem(1,2)*((-(r*((-t)+2*r-2))/2.0E+0)-(r*(1-t))/2.0E+0);
          Bt.xelem(2,1) = invJ.xelem(2,0)*(((1-t)*((-t)+2*r-2))/2.0E+0+r*(1-t))+invJ.xelem(2,2)*((-(r*((-t)+2*r-2))/2.0E+0)-(r*(1-t))/2.0E+0);
          Bt.xelem(0,2) = invJ.xelem(0,1)*(((1-t)*((-t)+2*s-2))/2.0E+0+s*(1-t))+invJ.xelem(0,2)*((-(s*((-t)+2*s-2))/2.0E+0)-(s*(1-t))/2.0E+0);
          Bt.xelem(1,2) = invJ.xelem(1,1)*(((1-t)*((-t)+2*s-2))/2.0E+0+s*(1-t))+invJ.xelem(1,2)*((-(s*((-t)+2*s-2))/2.0E+0)-(s*(1-t))/2.0E+0);
          Bt.xelem(2,2) = invJ.xelem(2,1)*(((1-t)*((-t)+2*s-2))/2.0E+0+s*(1-t))+invJ.xelem(2,2)*((-(s*((-t)+2*s-2))/2.0E+0)-(s*(1-t))/2.0E+0);
          Bt.xelem(0,3) = invJ.xelem(0,1)*((-((t+1)*(t-2*s-2*r))/2.0E+0)-((-s)-r+1)*(t+1))+invJ.xelem(0,0)*((-((t+1)*(t-2*s-2*r))/2.0E+0)-((-s)-r+1)*(t+1))+invJ.xelem(0,2)*((((-s)-r+1)*(t-2*s-2*r))/2.0E+0+(((-s)-r+1)*(t+1))/2.0E+0);
          Bt.xelem(1,3) = invJ.xelem(1,1)*((-((t+1)*(t-2*s-2*r))/2.0E+0)-((-s)-r+1)*(t+1))+invJ.xelem(1,0)*((-((t+1)*(t-2*s-2*r))/2.0E+0)-((-s)-r+1)*(t+1))+invJ.xelem(1,2)*((((-s)-r+1)*(t-2*s-2*r))/2.0E+0+(((-s)-r+1)*(t+1))/2.0E+0);
          Bt.xelem(2,3) = invJ.xelem(2,1)*((-((t+1)*(t-2*s-2*r))/2.0E+0)-((-s)-r+1)*(t+1))+invJ.xelem(2,0)*((-((t+1)*(t-2*s-2*r))/2.0E+0)-((-s)-r+1)*(t+1))+invJ.xelem(2,2)*((((-s)-r+1)*(t-2*s-2*r))/2.0E+0+(((-s)-r+1)*(t+1))/2.0E+0);
          Bt.xelem(0,4) = invJ.xelem(0,0)*(((t+1)*(t+2*r-2))/2.0E+0+r*(t+1))+invJ.xelem(0,2)*((r*(t+2*r-2))/2.0E+0+(r*(t+1))/2.0E+0);
          Bt.xelem(1,4) = invJ.xelem(1,0)*(((t+1)*(t+2*r-2))/2.0E+0+r*(t+1))+invJ.xelem(1,2)*((r*(t+2*r-2))/2.0E+0+(r*(t+1))/2.0E+0);
          Bt.xelem(2,4) = invJ.xelem(2,0)*(((t+1)*(t+2*r-2))/2.0E+0+r*(t+1))+invJ.xelem(2,2)*((r*(t+2*r-2))/2.0E+0+(r*(t+1))/2.0E+0);
          Bt.xelem(0,5) = invJ.xelem(0,1)*(((t+1)*(t+2*s-2))/2.0E+0+s*(t+1))+invJ.xelem(0,2)*((s*(t+2*s-2))/2.0E+0+(s*(t+1))/2.0E+0);
          Bt.xelem(1,5) = invJ.xelem(1,1)*(((t+1)*(t+2*s-2))/2.0E+0+s*(t+1))+invJ.xelem(1,2)*((s*(t+2*s-2))/2.0E+0+(s*(t+1))/2.0E+0);
          Bt.xelem(2,5) = invJ.xelem(2,1)*(((t+1)*(t+2*s-2))/2.0E+0+s*(t+1))+invJ.xelem(2,2)*((s*(t+2*s-2))/2.0E+0+(s*(t+1))/2.0E+0);
          Bt.xelem(0,6) = (-2*invJ.xelem(0,1)*r*(1-t))+invJ.xelem(0,0)*(2*((-s)-r+1)*(1-t)-2*r*(1-t))-2*invJ.xelem(0,2)*r*((-s)-r+1);
          Bt.xelem(1,6) = (-2*invJ.xelem(1,1)*r*(1-t))+invJ.xelem(1,0)*(2*((-s)-r+1)*(1-t)-2*r*(1-t))-2*invJ.xelem(1,2)*r*((-s)-r+1);
          Bt.xelem(2,6) = (-2*invJ.xelem(2,1)*r*(1-t))+invJ.xelem(2,0)*(2*((-s)-r+1)*(1-t)-2*r*(1-t))-2*invJ.xelem(2,2)*r*((-s)-r+1);
          Bt.xelem(0,7) = 2*invJ.xelem(0,0)*s*(1-t)+2*invJ.xelem(0,1)*r*(1-t)-2*invJ.xelem(0,2)*r*s;
          Bt.xelem(1,7) = 2*invJ.xelem(1,0)*s*(1-t)+2*invJ.xelem(1,1)*r*(1-t)-2*invJ.xelem(1,2)*r*s;
          Bt.xelem(2,7) = 2*invJ.xelem(2,0)*s*(1-t)+2*invJ.xelem(2,1)*r*(1-t)-2*invJ.xelem(2,2)*r*s;
          Bt.xelem(0,8) = (-2*invJ.xelem(0,0)*s*(1-t))+invJ.xelem(0,1)*(2*((-s)-r+1)*(1-t)-2*s*(1-t))-2*invJ.xelem(0,2)*((-s)-r+1)*s;
          Bt.xelem(1,8) = (-2*invJ.xelem(1,0)*s*(1-t))+invJ.xelem(1,1)*(2*((-s)-r+1)*(1-t)-2*s*(1-t))-2*invJ.xelem(1,2)*((-s)-r+1)*s;
          Bt.xelem(2,8) = (-2*invJ.xelem(2,0)*s*(1-t))+invJ.xelem(2,1)*(2*((-s)-r+1)*(1-t)-2*s*(1-t))-2*invJ.xelem(2,2)*((-s)-r+1)*s;
          Bt.xelem(0,9) = invJ.xelem(0,0)*(2*((-s)-r+1)*(t+1)-2*r*(t+1))-2*invJ.xelem(0,1)*r*(t+1)+2*invJ.xelem(0,2)*r*((-s)-r+1);
          Bt.xelem(1,9) = invJ.xelem(1,0)*(2*((-s)-r+1)*(t+1)-2*r*(t+1))-2*invJ.xelem(1,1)*r*(t+1)+2*invJ.xelem(1,2)*r*((-s)-r+1);
          Bt.xelem(2,9) = invJ.xelem(2,0)*(2*((-s)-r+1)*(t+1)-2*r*(t+1))-2*invJ.xelem(2,1)*r*(t+1)+2*invJ.xelem(2,2)*r*((-s)-r+1);
          Bt.xelem(0,10) = 2*invJ.xelem(0,0)*s*(t+1)+2*invJ.xelem(0,1)*r*(t+1)+2*invJ.xelem(0,2)*r*s;
          Bt.xelem(1,10) = 2*invJ.xelem(1,0)*s*(t+1)+2*invJ.xelem(1,1)*r*(t+1)+2*invJ.xelem(1,2)*r*s;
          Bt.xelem(2,10) = 2*invJ.xelem(2,0)*s*(t+1)+2*invJ.xelem(2,1)*r*(t+1)+2*invJ.xelem(2,2)*r*s;
          Bt.xelem(0,11) = invJ.xelem(0,1)*(2*((-s)-r+1)*(t+1)-2*s*(t+1))-2*invJ.xelem(0,0)*s*(t+1)+2*invJ.xelem(0,2)*((-s)-r+1)*s;
          Bt.xelem(1,11) = invJ.xelem(1,1)*(2*((-s)-r+1)*(t+1)-2*s*(t+1))-2*invJ.xelem(1,0)*s*(t+1)+2*invJ.xelem(1,2)*((-s)-r+1)*s;
          Bt.xelem(2,11) = invJ.xelem(2,1)*(2*((-s)-r+1)*(t+1)-2*s*(t+1))-2*invJ.xelem(2,0)*s*(t+1)+2*invJ.xelem(2,2)*((-s)-r+1)*s;
          Bt.xelem(0,12) = invJ.xelem(0,1)*(t_2-1)+invJ.xelem(0,0)*(t_2-1)-2*invJ.xelem(0,2)*((-s)-r+1)*t;
          Bt.xelem(1,12) = invJ.xelem(1,1)*(t_2-1)+invJ.xelem(1,0)*(t_2-1)-2*invJ.xelem(1,2)*((-s)-r+1)*t;
          Bt.xelem(2,12) = invJ.xelem(2,1)*(t_2-1)+invJ.xelem(2,0)*(t_2-1)-2*invJ.xelem(2,2)*((-s)-r+1)*t;
          Bt.xelem(0,13) = invJ.xelem(0,0)*(1-t_2)-2*invJ.xelem(0,2)*r*t;
          Bt.xelem(1,13) = invJ.xelem(1,0)*(1-t_2)-2*invJ.xelem(1,2)*r*t;
          Bt.xelem(2,13) = invJ.xelem(2,0)*(1-t_2)-2*invJ.xelem(2,2)*r*t;
          Bt.xelem(0,14) = invJ.xelem(0,1)*(1-t_2)-2*invJ.xelem(0,2)*s*t;
          Bt.xelem(1,14) = invJ.xelem(1,1)*(1-t_2)-2*invJ.xelem(1,2)*s*t;
          Bt.xelem(2,14) = invJ.xelem(2,1)*(1-t_2)-2*invJ.xelem(2,2)*s*t;
     }
     
     virtual Matrix InterpGaussToNodal(FemMatrixType eMatType, const Matrix& taug) const {
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

private:
     void ScalarInterpMatrix(const ColumnVector& rv, Matrix& Hs, octave_idx_type irow) const {
          FEM_ASSERT(rv.numel() == 3);
          FEM_ASSERT(Hs.columns() == 15);
          FEM_ASSERT(irow >= 0);
          FEM_ASSERT(irow < Hs.rows());

          const double r = rv.xelem(0);
          const double s = rv.xelem(1);
          const double t = rv.xelem(2);
          const double t2 = t * t;

          Hs.xelem(irow,0) = (((-s)-r+1)*(1-t)*((-t)-2*s-2*r))/2.0E+0;
          Hs.xelem(irow,1) = (r*(1-t)*((-t)+2*r-2))/2.0E+0;
          Hs.xelem(irow,2) = (s*(1-t)*((-t)+2*s-2))/2.0E+0;
          Hs.xelem(irow,3) = (((-s)-r+1)*(t+1)*(t-2*s-2*r))/2.0E+0;
          Hs.xelem(irow,4) = (r*(t+1)*(t+2*r-2))/2.0E+0;
          Hs.xelem(irow,5) = (s*(t+1)*(t+2*s-2))/2.0E+0;
          Hs.xelem(irow,6) = 2*r*((-s)-r+1)*(1-t);
          Hs.xelem(irow,7) = 2*r*s*(1-t);
          Hs.xelem(irow,8) = 2*((-s)-r+1)*s*(1-t);
          Hs.xelem(irow,9) = 2*r*((-s)-r+1)*(t+1);
          Hs.xelem(irow,10) = 2*r*s*(t+1);
          Hs.xelem(irow,11) = 2*((-s)-r+1)*s*(t+1);
          Hs.xelem(irow,12) = ((-s)-r+1)*(1-t2);
          Hs.xelem(irow,13) = r*(1-t2);
          Hs.xelem(irow,14) = s*(1-t2);
     }
};

class Tet10h: public Element3D
{
public:
     Tet10h(ElementTypes::TypeId eltype, octave_idx_type id, const Matrix& X, const Material* material, const int32NDArray& nodes, const StrainField& oRefStrain)
          :Element3D(eltype, id, X, material, nodes, oRefStrain) {
          FEM_ASSERT(nodes.numel() == 10);
     }

     virtual const IntegrationRule& GetIntegrationRule(FemMatrixType eMatType) const {

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

          enum IntegRuleType {
               R1 = 0,
               R2,
               R3,
               RNUM
          };

          static constexpr octave_idx_type Ni[RNUM] = {N1, N2, N3};
          static constexpr const double* ri[RNUM] = {r1, r2, r3};
          static constexpr const double* si[RNUM] = {s1, s2, s3};
          static constexpr const double* ti[RNUM] = {t1, t2, t3};
          static constexpr const double* wi[RNUM] = {w1, w2, w3};

          static array<IntegrationRule, RNUM> rgIntegRule;

          octave_idx_type iIntegRule;

          switch (eMatType) {
          case MAT_STIFFNESS:
          case MAT_STIFFNESS_SYM:
          case MAT_STIFFNESS_SYM_L:
          case VEC_LOAD_CONSISTENT:
          case VEC_LOAD_LUMPED:
          case MAT_THERMAL_COND:
               iIntegRule = R1;
               break;

          case MAT_MASS:
          case MAT_MASS_SYM:
          case MAT_MASS_SYM_L:
          case VEC_INERTIA_M1:
          case MAT_INERTIA_J:
          case MAT_INERTIA_INV3:
          case MAT_INERTIA_INV4:
          case MAT_INERTIA_INV5:
          case MAT_INERTIA_INV8:
          case MAT_INERTIA_INV9:
          case MAT_ACCEL_LOAD:
          case SCA_TOT_MASS:
          case MAT_HEAT_CAPACITY:
               iIntegRule = R2;
               break;

          case VEC_STRESS_CAUCH:
          case VEC_STRAIN_TOTAL:
          case SCA_STRESS_VMIS:
               iIntegRule = R3;
               break;

          default:
               throw std::runtime_error("matrix type not supported");
          }

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

          return rgIntegRule[iIntegRule];
     }

protected:
     virtual double Jacobian(const ColumnVector& rv, Matrix& J) const {
          FEM_ASSERT(J.rows() == 3);
          FEM_ASSERT(J.columns() == 3);
          FEM_ASSERT(rv.numel() == 3);

          const double r = rv.xelem(0);
          const double s = rv.xelem(1);
          const double t = rv.xelem(2);

          J.xelem(0,0) = 4*X.xelem(0,8)*t-4*X.xelem(0,5)*t+4*X.xelem(0,9)*((-t)-s-r+1)-2*X.xelem(0,2)*((-t)-s-r+1)-X.xelem(0,2)*((-2*t)-2*s-2*r+1)+4*X.xelem(0,7)*s-4*X.xelem(0,6)*s+X.xelem(0,3)*(2*r-1)-4*X.xelem(0,9)*r+2*X.xelem(0,3)*r;
          J.xelem(1,0) = (-4*X.xelem(0,5)*t)+4*X.xelem(0,4)*t+4*X.xelem(0,6)*((-t)-s-r+1)-2*X.xelem(0,2)*((-t)-s-r+1)-X.xelem(0,2)*((-2*t)-2*s-2*r+1)+X.xelem(0,0)*(2*s-1)-4*X.xelem(0,6)*s+2*X.xelem(0,0)*s-4*X.xelem(0,9)*r+4*X.xelem(0,7)*r;
          J.xelem(2,0) = X.xelem(0,1)*(2*t-1)-4*X.xelem(0,5)*t+2*X.xelem(0,1)*t+4*X.xelem(0,5)*((-t)-s-r+1)-2*X.xelem(0,2)*((-t)-s-r+1)-X.xelem(0,2)*((-2*t)-2*s-2*r+1)-4*X.xelem(0,6)*s+4*X.xelem(0,4)*s-4*X.xelem(0,9)*r+4*X.xelem(0,8)*r;
          J.xelem(0,1) = 4*X.xelem(1,8)*t-4*X.xelem(1,5)*t+4*X.xelem(1,9)*((-t)-s-r+1)-2*X.xelem(1,2)*((-t)-s-r+1)-X.xelem(1,2)*((-2*t)-2*s-2*r+1)+4*X.xelem(1,7)*s-4*X.xelem(1,6)*s+X.xelem(1,3)*(2*r-1)-4*X.xelem(1,9)*r+2*X.xelem(1,3)*r;
          J.xelem(1,1) = (-4*X.xelem(1,5)*t)+4*X.xelem(1,4)*t+4*X.xelem(1,6)*((-t)-s-r+1)-2*X.xelem(1,2)*((-t)-s-r+1)-X.xelem(1,2)*((-2*t)-2*s-2*r+1)+X.xelem(1,0)*(2*s-1)-4*X.xelem(1,6)*s+2*X.xelem(1,0)*s-4*X.xelem(1,9)*r+4*X.xelem(1,7)*r;
          J.xelem(2,1) = X.xelem(1,1)*(2*t-1)-4*X.xelem(1,5)*t+2*X.xelem(1,1)*t+4*X.xelem(1,5)*((-t)-s-r+1)-2*X.xelem(1,2)*((-t)-s-r+1)-X.xelem(1,2)*((-2*t)-2*s-2*r+1)-4*X.xelem(1,6)*s+4*X.xelem(1,4)*s-4*X.xelem(1,9)*r+4*X.xelem(1,8)*r;
          J.xelem(0,2) = 4*X.xelem(2,8)*t-4*X.xelem(2,5)*t+4*X.xelem(2,9)*((-t)-s-r+1)-2*X.xelem(2,2)*((-t)-s-r+1)-X.xelem(2,2)*((-2*t)-2*s-2*r+1)+4*X.xelem(2,7)*s-4*X.xelem(2,6)*s+X.xelem(2,3)*(2*r-1)-4*X.xelem(2,9)*r+2*X.xelem(2,3)*r;
          J.xelem(1,2) = (-4*X.xelem(2,5)*t)+4*X.xelem(2,4)*t+4*X.xelem(2,6)*((-t)-s-r+1)-2*X.xelem(2,2)*((-t)-s-r+1)-X.xelem(2,2)*((-2*t)-2*s-2*r+1)+X.xelem(2,0)*(2*s-1)-4*X.xelem(2,6)*s+2*X.xelem(2,0)*s-4*X.xelem(2,9)*r+4*X.xelem(2,7)*r;
          J.xelem(2,2) = X.xelem(2,1)*(2*t-1)-4*X.xelem(2,5)*t+2*X.xelem(2,1)*t+4*X.xelem(2,5)*((-t)-s-r+1)-2*X.xelem(2,2)*((-t)-s-r+1)-X.xelem(2,2)*((-2*t)-2*s-2*r+1)-4*X.xelem(2,6)*s+4*X.xelem(2,4)*s-4*X.xelem(2,9)*r+4*X.xelem(2,8)*r;

          return J.determinant();
     }

     virtual void DispInterpMatrix(const ColumnVector& rv, Matrix& H) const {
          FEM_ASSERT(rv.numel() == 3);
          FEM_ASSERT(H.rows() == 3);
          FEM_ASSERT(H.columns() == 30);

          const double r = rv.xelem(0);
          const double s = rv.xelem(1);
          const double t = rv.xelem(2);

          H.xelem(0,0) = s*(2*s-1);
          H.xelem(1,0) = 0;
          H.xelem(2,0) = 0;
          H.xelem(0,1) = 0;
          H.xelem(1,1) = s*(2*s-1);
          H.xelem(2,1) = 0;
          H.xelem(0,2) = 0;
          H.xelem(1,2) = 0;
          H.xelem(2,2) = s*(2*s-1);
          H.xelem(0,3) = t*(2*t-1);
          H.xelem(1,3) = 0;
          H.xelem(2,3) = 0;
          H.xelem(0,4) = 0;
          H.xelem(1,4) = t*(2*t-1);
          H.xelem(2,4) = 0;
          H.xelem(0,5) = 0;
          H.xelem(1,5) = 0;
          H.xelem(2,5) = t*(2*t-1);
          H.xelem(0,6) = ((-2*t)-2*s-2*r+1)*((-t)-s-r+1);
          H.xelem(1,6) = 0;
          H.xelem(2,6) = 0;
          H.xelem(0,7) = 0;
          H.xelem(1,7) = ((-2*t)-2*s-2*r+1)*((-t)-s-r+1);
          H.xelem(2,7) = 0;
          H.xelem(0,8) = 0;
          H.xelem(1,8) = 0;
          H.xelem(2,8) = ((-2*t)-2*s-2*r+1)*((-t)-s-r+1);
          H.xelem(0,9) = r*(2*r-1);
          H.xelem(1,9) = 0;
          H.xelem(2,9) = 0;
          H.xelem(0,10) = 0;
          H.xelem(1,10) = r*(2*r-1);
          H.xelem(2,10) = 0;
          H.xelem(0,11) = 0;
          H.xelem(1,11) = 0;
          H.xelem(2,11) = r*(2*r-1);
          H.xelem(0,12) = 4*s*t;
          H.xelem(1,12) = 0;
          H.xelem(2,12) = 0;
          H.xelem(0,13) = 0;
          H.xelem(1,13) = 4*s*t;
          H.xelem(2,13) = 0;
          H.xelem(0,14) = 0;
          H.xelem(1,14) = 0;
          H.xelem(2,14) = 4*s*t;
          H.xelem(0,15) = 4*((-t)-s-r+1)*t;
          H.xelem(1,15) = 0;
          H.xelem(2,15) = 0;
          H.xelem(0,16) = 0;
          H.xelem(1,16) = 4*((-t)-s-r+1)*t;
          H.xelem(2,16) = 0;
          H.xelem(0,17) = 0;
          H.xelem(1,17) = 0;
          H.xelem(2,17) = 4*((-t)-s-r+1)*t;
          H.xelem(0,18) = 4*s*((-t)-s-r+1);
          H.xelem(1,18) = 0;
          H.xelem(2,18) = 0;
          H.xelem(0,19) = 0;
          H.xelem(1,19) = 4*s*((-t)-s-r+1);
          H.xelem(2,19) = 0;
          H.xelem(0,20) = 0;
          H.xelem(1,20) = 0;
          H.xelem(2,20) = 4*s*((-t)-s-r+1);
          H.xelem(0,21) = 4*r*s;
          H.xelem(1,21) = 0;
          H.xelem(2,21) = 0;
          H.xelem(0,22) = 0;
          H.xelem(1,22) = 4*r*s;
          H.xelem(2,22) = 0;
          H.xelem(0,23) = 0;
          H.xelem(1,23) = 0;
          H.xelem(2,23) = 4*r*s;
          H.xelem(0,24) = 4*r*t;
          H.xelem(1,24) = 0;
          H.xelem(2,24) = 0;
          H.xelem(0,25) = 0;
          H.xelem(1,25) = 4*r*t;
          H.xelem(2,25) = 0;
          H.xelem(0,26) = 0;
          H.xelem(1,26) = 0;
          H.xelem(2,26) = 4*r*t;
          H.xelem(0,27) = 4*r*((-t)-s-r+1);
          H.xelem(1,27) = 0;
          H.xelem(2,27) = 0;
          H.xelem(0,28) = 0;
          H.xelem(1,28) = 4*r*((-t)-s-r+1);
          H.xelem(2,28) = 0;
          H.xelem(0,29) = 0;
          H.xelem(1,29) = 0;
          H.xelem(2,29) = 4*r*((-t)-s-r+1);
     }

     virtual void StrainMatrix(const ColumnVector& rv, const Matrix& invJ, Matrix& B) const {
          FEM_ASSERT(rv.numel() == 3);
          FEM_ASSERT(invJ.rows() == 3);
          FEM_ASSERT(invJ.columns() == 3);
          FEM_ASSERT(B.rows() == 6);
          FEM_ASSERT(B.columns() == 30);

          const double r = rv.xelem(0);
          const double s = rv.xelem(1);
          const double t = rv.xelem(2);

          B.xelem(0,0) = invJ.xelem(0,1)*(4*s-1);
          B.xelem(1,0) = 0;
          B.xelem(2,0) = 0;
          B.xelem(3,0) = invJ.xelem(1,1)*(4*s-1);
          B.xelem(4,0) = 0;
          B.xelem(5,0) = invJ.xelem(2,1)*(4*s-1);
          B.xelem(0,1) = 0;
          B.xelem(1,1) = invJ.xelem(1,1)*(4*s-1);
          B.xelem(2,1) = 0;
          B.xelem(3,1) = invJ.xelem(0,1)*(4*s-1);
          B.xelem(4,1) = invJ.xelem(2,1)*(4*s-1);
          B.xelem(5,1) = 0;
          B.xelem(0,2) = 0;
          B.xelem(1,2) = 0;
          B.xelem(2,2) = invJ.xelem(2,1)*(4*s-1);
          B.xelem(3,2) = 0;
          B.xelem(4,2) = invJ.xelem(1,1)*(4*s-1);
          B.xelem(5,2) = invJ.xelem(0,1)*(4*s-1);
          B.xelem(0,3) = invJ.xelem(0,2)*(4*t-1);
          B.xelem(1,3) = 0;
          B.xelem(2,3) = 0;
          B.xelem(3,3) = invJ.xelem(1,2)*(4*t-1);
          B.xelem(4,3) = 0;
          B.xelem(5,3) = invJ.xelem(2,2)*(4*t-1);
          B.xelem(0,4) = 0;
          B.xelem(1,4) = invJ.xelem(1,2)*(4*t-1);
          B.xelem(2,4) = 0;
          B.xelem(3,4) = invJ.xelem(0,2)*(4*t-1);
          B.xelem(4,4) = invJ.xelem(2,2)*(4*t-1);
          B.xelem(5,4) = 0;
          B.xelem(0,5) = 0;
          B.xelem(1,5) = 0;
          B.xelem(2,5) = invJ.xelem(2,2)*(4*t-1);
          B.xelem(3,5) = 0;
          B.xelem(4,5) = invJ.xelem(1,2)*(4*t-1);
          B.xelem(5,5) = invJ.xelem(0,2)*(4*t-1);
          B.xelem(0,6) = invJ.xelem(0,2)*(2*t-2*((-t)-s-r+1)+2*s+2*r-1)+invJ.xelem(0,1)*(2*t-2*((-t)-s-r+1)+2*s+2*r-1)+invJ.xelem(0,0)*(2*t-2*((-t)-s-r+1)+2*s+2*r-1);
          B.xelem(1,6) = 0;
          B.xelem(2,6) = 0;
          B.xelem(3,6) = invJ.xelem(1,2)*(2*t-2*((-t)-s-r+1)+2*s+2*r-1)+invJ.xelem(1,1)*(2*t-2*((-t)-s-r+1)+2*s+2*r-1)+invJ.xelem(1,0)*(2*t-2*((-t)-s-r+1)+2*s+2*r-1);
          B.xelem(4,6) = 0;
          B.xelem(5,6) = invJ.xelem(2,2)*(2*t-2*((-t)-s-r+1)+2*s+2*r-1)+invJ.xelem(2,1)*(2*t-2*((-t)-s-r+1)+2*s+2*r-1)+invJ.xelem(2,0)*(2*t-2*((-t)-s-r+1)+2*s+2*r-1);
          B.xelem(0,7) = 0;
          B.xelem(1,7) = invJ.xelem(1,2)*(2*t-2*((-t)-s-r+1)+2*s+2*r-1)+invJ.xelem(1,1)*(2*t-2*((-t)-s-r+1)+2*s+2*r-1)+invJ.xelem(1,0)*(2*t-2*((-t)-s-r+1)+2*s+2*r-1);
          B.xelem(2,7) = 0;
          B.xelem(3,7) = invJ.xelem(0,2)*(2*t-2*((-t)-s-r+1)+2*s+2*r-1)+invJ.xelem(0,1)*(2*t-2*((-t)-s-r+1)+2*s+2*r-1)+invJ.xelem(0,0)*(2*t-2*((-t)-s-r+1)+2*s+2*r-1);
          B.xelem(4,7) = invJ.xelem(2,2)*(2*t-2*((-t)-s-r+1)+2*s+2*r-1)+invJ.xelem(2,1)*(2*t-2*((-t)-s-r+1)+2*s+2*r-1)+invJ.xelem(2,0)*(2*t-2*((-t)-s-r+1)+2*s+2*r-1);
          B.xelem(5,7) = 0;
          B.xelem(0,8) = 0;
          B.xelem(1,8) = 0;
          B.xelem(2,8) = invJ.xelem(2,2)*(2*t-2*((-t)-s-r+1)+2*s+2*r-1)+invJ.xelem(2,1)*(2*t-2*((-t)-s-r+1)+2*s+2*r-1)+invJ.xelem(2,0)*(2*t-2*((-t)-s-r+1)+2*s+2*r-1);
          B.xelem(3,8) = 0;
          B.xelem(4,8) = invJ.xelem(1,2)*(2*t-2*((-t)-s-r+1)+2*s+2*r-1)+invJ.xelem(1,1)*(2*t-2*((-t)-s-r+1)+2*s+2*r-1)+invJ.xelem(1,0)*(2*t-2*((-t)-s-r+1)+2*s+2*r-1);
          B.xelem(5,8) = invJ.xelem(0,2)*(2*t-2*((-t)-s-r+1)+2*s+2*r-1)+invJ.xelem(0,1)*(2*t-2*((-t)-s-r+1)+2*s+2*r-1)+invJ.xelem(0,0)*(2*t-2*((-t)-s-r+1)+2*s+2*r-1);
          B.xelem(0,9) = invJ.xelem(0,0)*(4*r-1);
          B.xelem(1,9) = 0;
          B.xelem(2,9) = 0;
          B.xelem(3,9) = invJ.xelem(1,0)*(4*r-1);
          B.xelem(4,9) = 0;
          B.xelem(5,9) = invJ.xelem(2,0)*(4*r-1);
          B.xelem(0,10) = 0;
          B.xelem(1,10) = invJ.xelem(1,0)*(4*r-1);
          B.xelem(2,10) = 0;
          B.xelem(3,10) = invJ.xelem(0,0)*(4*r-1);
          B.xelem(4,10) = invJ.xelem(2,0)*(4*r-1);
          B.xelem(5,10) = 0;
          B.xelem(0,11) = 0;
          B.xelem(1,11) = 0;
          B.xelem(2,11) = invJ.xelem(2,0)*(4*r-1);
          B.xelem(3,11) = 0;
          B.xelem(4,11) = invJ.xelem(1,0)*(4*r-1);
          B.xelem(5,11) = invJ.xelem(0,0)*(4*r-1);
          B.xelem(0,12) = 4*invJ.xelem(0,1)*t+4*invJ.xelem(0,2)*s;
          B.xelem(1,12) = 0;
          B.xelem(2,12) = 0;
          B.xelem(3,12) = 4*invJ.xelem(1,1)*t+4*invJ.xelem(1,2)*s;
          B.xelem(4,12) = 0;
          B.xelem(5,12) = 4*invJ.xelem(2,1)*t+4*invJ.xelem(2,2)*s;
          B.xelem(0,13) = 0;
          B.xelem(1,13) = 4*invJ.xelem(1,1)*t+4*invJ.xelem(1,2)*s;
          B.xelem(2,13) = 0;
          B.xelem(3,13) = 4*invJ.xelem(0,1)*t+4*invJ.xelem(0,2)*s;
          B.xelem(4,13) = 4*invJ.xelem(2,1)*t+4*invJ.xelem(2,2)*s;
          B.xelem(5,13) = 0;
          B.xelem(0,14) = 0;
          B.xelem(1,14) = 0;
          B.xelem(2,14) = 4*invJ.xelem(2,1)*t+4*invJ.xelem(2,2)*s;
          B.xelem(3,14) = 0;
          B.xelem(4,14) = 4*invJ.xelem(1,1)*t+4*invJ.xelem(1,2)*s;
          B.xelem(5,14) = 4*invJ.xelem(0,1)*t+4*invJ.xelem(0,2)*s;
          B.xelem(0,15) = (-4*invJ.xelem(0,1)*t)-4*invJ.xelem(0,0)*t+invJ.xelem(0,2)*(4*((-t)-s-r+1)-4*t);
          B.xelem(1,15) = 0;
          B.xelem(2,15) = 0;
          B.xelem(3,15) = (-4*invJ.xelem(1,1)*t)-4*invJ.xelem(1,0)*t+invJ.xelem(1,2)*(4*((-t)-s-r+1)-4*t);
          B.xelem(4,15) = 0;
          B.xelem(5,15) = (-4*invJ.xelem(2,1)*t)-4*invJ.xelem(2,0)*t+invJ.xelem(2,2)*(4*((-t)-s-r+1)-4*t);
          B.xelem(0,16) = 0;
          B.xelem(1,16) = (-4*invJ.xelem(1,1)*t)-4*invJ.xelem(1,0)*t+invJ.xelem(1,2)*(4*((-t)-s-r+1)-4*t);
          B.xelem(2,16) = 0;
          B.xelem(3,16) = (-4*invJ.xelem(0,1)*t)-4*invJ.xelem(0,0)*t+invJ.xelem(0,2)*(4*((-t)-s-r+1)-4*t);
          B.xelem(4,16) = (-4*invJ.xelem(2,1)*t)-4*invJ.xelem(2,0)*t+invJ.xelem(2,2)*(4*((-t)-s-r+1)-4*t);
          B.xelem(5,16) = 0;
          B.xelem(0,17) = 0;
          B.xelem(1,17) = 0;
          B.xelem(2,17) = (-4*invJ.xelem(2,1)*t)-4*invJ.xelem(2,0)*t+invJ.xelem(2,2)*(4*((-t)-s-r+1)-4*t);
          B.xelem(3,17) = 0;
          B.xelem(4,17) = (-4*invJ.xelem(1,1)*t)-4*invJ.xelem(1,0)*t+invJ.xelem(1,2)*(4*((-t)-s-r+1)-4*t);
          B.xelem(5,17) = (-4*invJ.xelem(0,1)*t)-4*invJ.xelem(0,0)*t+invJ.xelem(0,2)*(4*((-t)-s-r+1)-4*t);
          B.xelem(0,18) = invJ.xelem(0,1)*(4*((-t)-s-r+1)-4*s)-4*invJ.xelem(0,2)*s-4*invJ.xelem(0,0)*s;
          B.xelem(1,18) = 0;
          B.xelem(2,18) = 0;
          B.xelem(3,18) = invJ.xelem(1,1)*(4*((-t)-s-r+1)-4*s)-4*invJ.xelem(1,2)*s-4*invJ.xelem(1,0)*s;
          B.xelem(4,18) = 0;
          B.xelem(5,18) = invJ.xelem(2,1)*(4*((-t)-s-r+1)-4*s)-4*invJ.xelem(2,2)*s-4*invJ.xelem(2,0)*s;
          B.xelem(0,19) = 0;
          B.xelem(1,19) = invJ.xelem(1,1)*(4*((-t)-s-r+1)-4*s)-4*invJ.xelem(1,2)*s-4*invJ.xelem(1,0)*s;
          B.xelem(2,19) = 0;
          B.xelem(3,19) = invJ.xelem(0,1)*(4*((-t)-s-r+1)-4*s)-4*invJ.xelem(0,2)*s-4*invJ.xelem(0,0)*s;
          B.xelem(4,19) = invJ.xelem(2,1)*(4*((-t)-s-r+1)-4*s)-4*invJ.xelem(2,2)*s-4*invJ.xelem(2,0)*s;
          B.xelem(5,19) = 0;
          B.xelem(0,20) = 0;
          B.xelem(1,20) = 0;
          B.xelem(2,20) = invJ.xelem(2,1)*(4*((-t)-s-r+1)-4*s)-4*invJ.xelem(2,2)*s-4*invJ.xelem(2,0)*s;
          B.xelem(3,20) = 0;
          B.xelem(4,20) = invJ.xelem(1,1)*(4*((-t)-s-r+1)-4*s)-4*invJ.xelem(1,2)*s-4*invJ.xelem(1,0)*s;
          B.xelem(5,20) = invJ.xelem(0,1)*(4*((-t)-s-r+1)-4*s)-4*invJ.xelem(0,2)*s-4*invJ.xelem(0,0)*s;
          B.xelem(0,21) = 4*invJ.xelem(0,0)*s+4*invJ.xelem(0,1)*r;
          B.xelem(1,21) = 0;
          B.xelem(2,21) = 0;
          B.xelem(3,21) = 4*invJ.xelem(1,0)*s+4*invJ.xelem(1,1)*r;
          B.xelem(4,21) = 0;
          B.xelem(5,21) = 4*invJ.xelem(2,0)*s+4*invJ.xelem(2,1)*r;
          B.xelem(0,22) = 0;
          B.xelem(1,22) = 4*invJ.xelem(1,0)*s+4*invJ.xelem(1,1)*r;
          B.xelem(2,22) = 0;
          B.xelem(3,22) = 4*invJ.xelem(0,0)*s+4*invJ.xelem(0,1)*r;
          B.xelem(4,22) = 4*invJ.xelem(2,0)*s+4*invJ.xelem(2,1)*r;
          B.xelem(5,22) = 0;
          B.xelem(0,23) = 0;
          B.xelem(1,23) = 0;
          B.xelem(2,23) = 4*invJ.xelem(2,0)*s+4*invJ.xelem(2,1)*r;
          B.xelem(3,23) = 0;
          B.xelem(4,23) = 4*invJ.xelem(1,0)*s+4*invJ.xelem(1,1)*r;
          B.xelem(5,23) = 4*invJ.xelem(0,0)*s+4*invJ.xelem(0,1)*r;
          B.xelem(0,24) = 4*invJ.xelem(0,0)*t+4*invJ.xelem(0,2)*r;
          B.xelem(1,24) = 0;
          B.xelem(2,24) = 0;
          B.xelem(3,24) = 4*invJ.xelem(1,0)*t+4*invJ.xelem(1,2)*r;
          B.xelem(4,24) = 0;
          B.xelem(5,24) = 4*invJ.xelem(2,0)*t+4*invJ.xelem(2,2)*r;
          B.xelem(0,25) = 0;
          B.xelem(1,25) = 4*invJ.xelem(1,0)*t+4*invJ.xelem(1,2)*r;
          B.xelem(2,25) = 0;
          B.xelem(3,25) = 4*invJ.xelem(0,0)*t+4*invJ.xelem(0,2)*r;
          B.xelem(4,25) = 4*invJ.xelem(2,0)*t+4*invJ.xelem(2,2)*r;
          B.xelem(5,25) = 0;
          B.xelem(0,26) = 0;
          B.xelem(1,26) = 0;
          B.xelem(2,26) = 4*invJ.xelem(2,0)*t+4*invJ.xelem(2,2)*r;
          B.xelem(3,26) = 0;
          B.xelem(4,26) = 4*invJ.xelem(1,0)*t+4*invJ.xelem(1,2)*r;
          B.xelem(5,26) = 4*invJ.xelem(0,0)*t+4*invJ.xelem(0,2)*r;
          B.xelem(0,27) = invJ.xelem(0,0)*(4*((-t)-s-r+1)-4*r)-4*invJ.xelem(0,2)*r-4*invJ.xelem(0,1)*r;
          B.xelem(1,27) = 0;
          B.xelem(2,27) = 0;
          B.xelem(3,27) = invJ.xelem(1,0)*(4*((-t)-s-r+1)-4*r)-4*invJ.xelem(1,2)*r-4*invJ.xelem(1,1)*r;
          B.xelem(4,27) = 0;
          B.xelem(5,27) = invJ.xelem(2,0)*(4*((-t)-s-r+1)-4*r)-4*invJ.xelem(2,2)*r-4*invJ.xelem(2,1)*r;
          B.xelem(0,28) = 0;
          B.xelem(1,28) = invJ.xelem(1,0)*(4*((-t)-s-r+1)-4*r)-4*invJ.xelem(1,2)*r-4*invJ.xelem(1,1)*r;
          B.xelem(2,28) = 0;
          B.xelem(3,28) = invJ.xelem(0,0)*(4*((-t)-s-r+1)-4*r)-4*invJ.xelem(0,2)*r-4*invJ.xelem(0,1)*r;
          B.xelem(4,28) = invJ.xelem(2,0)*(4*((-t)-s-r+1)-4*r)-4*invJ.xelem(2,2)*r-4*invJ.xelem(2,1)*r;
          B.xelem(5,28) = 0;
          B.xelem(0,29) = 0;
          B.xelem(1,29) = 0;
          B.xelem(2,29) = invJ.xelem(2,0)*(4*((-t)-s-r+1)-4*r)-4*invJ.xelem(2,2)*r-4*invJ.xelem(2,1)*r;
          B.xelem(3,29) = 0;
          B.xelem(4,29) = invJ.xelem(1,0)*(4*((-t)-s-r+1)-4*r)-4*invJ.xelem(1,2)*r-4*invJ.xelem(1,1)*r;
          B.xelem(5,29) = invJ.xelem(0,0)*(4*((-t)-s-r+1)-4*r)-4*invJ.xelem(0,2)*r-4*invJ.xelem(0,1)*r;
     }

     virtual Matrix InterpGaussToNodal(FemMatrixType eMatType, const Matrix& taug) const {
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

          Matrix taun = H.solve(taug);

          taun.resize(nodes.numel(), taun.columns());

          static const struct {
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

     void ScalarInterpMatrix(const ColumnVector& rv, Matrix& Hs, octave_idx_type irow) const {
          FEM_ASSERT(rv.numel() == 3);
          FEM_ASSERT(Hs.columns() == 10);
          FEM_ASSERT(irow >= 0);
          FEM_ASSERT(irow < Hs.rows());

          const double r = rv.xelem(0);
          const double s = rv.xelem(1);
          const double t = rv.xelem(2);
          
          Hs.xelem(irow,0) = s*(2*s-1);
          Hs.xelem(irow,1) = t*(2*t-1);
          Hs.xelem(irow,2) = ((-2*t)-2*s-2*r+1)*((-t)-s-r+1);
          Hs.xelem(irow,3) = r*(2*r-1);
          Hs.xelem(irow,4) = 4*s*t;
          Hs.xelem(irow,5) = 4*((-t)-s-r+1)*t;
          Hs.xelem(irow,6) = 4*s*((-t)-s-r+1);
          Hs.xelem(irow,7) = 4*r*s;
          Hs.xelem(irow,8) = 4*r*t;
          Hs.xelem(irow,9) = 4*r*((-t)-s-r+1);
     }

     virtual void TemperatureGradientMatrix(const ColumnVector& rv, const Matrix& invJ, Matrix& Bt) const {
          FEM_ASSERT(Bt.rows() == 3);
          FEM_ASSERT(Bt.columns() == 10);
          FEM_ASSERT(rv.numel() == 3);

          const double r = rv.xelem(0);
          const double s = rv.xelem(1);
          const double t = rv.xelem(2);

          Bt.xelem(0,0) = invJ.xelem(0,1)*(4*s-1);
          Bt.xelem(1,0) = invJ.xelem(1,1)*(4*s-1);
          Bt.xelem(2,0) = invJ.xelem(2,1)*(4*s-1);
          Bt.xelem(0,1) = invJ.xelem(0,2)*(4*t-1);
          Bt.xelem(1,1) = invJ.xelem(1,2)*(4*t-1);
          Bt.xelem(2,1) = invJ.xelem(2,2)*(4*t-1);
          Bt.xelem(0,2) = invJ.xelem(0,2)*(2*t-2*((-t)-s-r+1)+2*s+2*r-1)+invJ.xelem(0,1)*(2*t-2*((-t)-s-r+1)+2*s+2*r-1)+invJ.xelem(0,0)*(2*t-2*((-t)-s-r+1)+2*s+2*r-1);
          Bt.xelem(1,2) = invJ.xelem(1,2)*(2*t-2*((-t)-s-r+1)+2*s+2*r-1)+invJ.xelem(1,1)*(2*t-2*((-t)-s-r+1)+2*s+2*r-1)+invJ.xelem(1,0)*(2*t-2*((-t)-s-r+1)+2*s+2*r-1);
          Bt.xelem(2,2) = invJ.xelem(2,2)*(2*t-2*((-t)-s-r+1)+2*s+2*r-1)+invJ.xelem(2,1)*(2*t-2*((-t)-s-r+1)+2*s+2*r-1)+invJ.xelem(2,0)*(2*t-2*((-t)-s-r+1)+2*s+2*r-1);
          Bt.xelem(0,3) = invJ.xelem(0,0)*(4*r-1);
          Bt.xelem(1,3) = invJ.xelem(1,0)*(4*r-1);
          Bt.xelem(2,3) = invJ.xelem(2,0)*(4*r-1);
          Bt.xelem(0,4) = 4*invJ.xelem(0,1)*t+4*invJ.xelem(0,2)*s;
          Bt.xelem(1,4) = 4*invJ.xelem(1,1)*t+4*invJ.xelem(1,2)*s;
          Bt.xelem(2,4) = 4*invJ.xelem(2,1)*t+4*invJ.xelem(2,2)*s;
          Bt.xelem(0,5) = (-4*invJ.xelem(0,1)*t)-4*invJ.xelem(0,0)*t+invJ.xelem(0,2)*(4*((-t)-s-r+1)-4*t);
          Bt.xelem(1,5) = (-4*invJ.xelem(1,1)*t)-4*invJ.xelem(1,0)*t+invJ.xelem(1,2)*(4*((-t)-s-r+1)-4*t);
          Bt.xelem(2,5) = (-4*invJ.xelem(2,1)*t)-4*invJ.xelem(2,0)*t+invJ.xelem(2,2)*(4*((-t)-s-r+1)-4*t);
          Bt.xelem(0,6) = invJ.xelem(0,1)*(4*((-t)-s-r+1)-4*s)-4*invJ.xelem(0,2)*s-4*invJ.xelem(0,0)*s;
          Bt.xelem(1,6) = invJ.xelem(1,1)*(4*((-t)-s-r+1)-4*s)-4*invJ.xelem(1,2)*s-4*invJ.xelem(1,0)*s;
          Bt.xelem(2,6) = invJ.xelem(2,1)*(4*((-t)-s-r+1)-4*s)-4*invJ.xelem(2,2)*s-4*invJ.xelem(2,0)*s;
          Bt.xelem(0,7) = 4*invJ.xelem(0,0)*s+4*invJ.xelem(0,1)*r;
          Bt.xelem(1,7) = 4*invJ.xelem(1,0)*s+4*invJ.xelem(1,1)*r;
          Bt.xelem(2,7) = 4*invJ.xelem(2,0)*s+4*invJ.xelem(2,1)*r;
          Bt.xelem(0,8) = 4*invJ.xelem(0,0)*t+4*invJ.xelem(0,2)*r;
          Bt.xelem(1,8) = 4*invJ.xelem(1,0)*t+4*invJ.xelem(1,2)*r;
          Bt.xelem(2,8) = 4*invJ.xelem(2,0)*t+4*invJ.xelem(2,2)*r;
          Bt.xelem(0,9) = invJ.xelem(0,0)*(4*((-t)-s-r+1)-4*r)-4*invJ.xelem(0,2)*r-4*invJ.xelem(0,1)*r;
          Bt.xelem(1,9) = invJ.xelem(1,0)*(4*((-t)-s-r+1)-4*r)-4*invJ.xelem(1,2)*r-4*invJ.xelem(1,1)*r;
          Bt.xelem(2,9) = invJ.xelem(2,0)*(4*((-t)-s-r+1)-4*r)-4*invJ.xelem(2,2)*r-4*invJ.xelem(2,1)*r;
     }
     
private:
     void ScalarInterpMatrixCornerNodes(const ColumnVector& rv, Matrix& Hs, octave_idx_type irow) const {
          FEM_ASSERT(rv.numel() == 3);
          FEM_ASSERT(Hs.columns() == 4);
          FEM_ASSERT(irow >= 0);
          FEM_ASSERT(irow < Hs.rows());

          const double r = rv.xelem(0);
          const double s = rv.xelem(1);
          const double t = rv.xelem(2);

          Hs.xelem(irow, 0) = s;
          Hs.xelem(irow, 1) = t;
          Hs.xelem(irow, 2) = 1. - r - s - t;
          Hs.xelem(irow, 3) = r;
     }
};

class Tet10: public Element3D
{
     static constexpr double gamma = 1. / 6.;

public:
     Tet10(ElementTypes::TypeId eltype, octave_idx_type id, const Matrix& X, const Material* material, const int32NDArray& nodes, const StrainField& oRefStrain)
          :Element3D(eltype, id, X, material, nodes, oRefStrain) {
          FEM_ASSERT(nodes.numel() == 10);
     }

     virtual const IntegrationRule& GetIntegrationRule(FemMatrixType eMatType) const {
          static IntegrationRule oIntegStiff, oIntegMass, oIntegMassDiag;

          switch (eMatType) {
          case MAT_STIFFNESS:
          case MAT_STIFFNESS_SYM:
          case MAT_STIFFNESS_SYM_L:
          case VEC_STRESS_CAUCH:
          case VEC_STRAIN_TOTAL:
          case VEC_LOAD_CONSISTENT:
          case VEC_LOAD_LUMPED:
          case MAT_THERMAL_COND:
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

               return oIntegStiff;

          case MAT_MASS:
          case MAT_MASS_SYM:
          case MAT_MASS_SYM_L:
          case VEC_INERTIA_M1:
          case MAT_INERTIA_J:
          case MAT_INERTIA_INV3:
          case MAT_INERTIA_INV4:
          case MAT_INERTIA_INV5:
          case MAT_INERTIA_INV8:
          case MAT_INERTIA_INV9:
          case MAT_ACCEL_LOAD:
          case MAT_HEAT_CAPACITY:
               if (!oIntegMass.iGetNumEvalPoints()) {
                    constexpr double g1 = 0.09273525031089122640232391373703060;
                    constexpr double g2 = 0.31088591926330060979734573376345783;
                    constexpr double g3 = 0.45449629587435035050811947372066056;
                    constexpr double w1 = (-1+6*g2*(2+g2*(-7+8*g2))+14*g3-60*g2*(3+4*g2*(-3+4*g2))*g3+4*(-7+30*g2*(3+4*g2*(-3+4*g2)))*g3*g3)/(120*(g1-g2)*(g2*(-3+8*g2)+6*g3+8*g2*(-3+4*g2)*g3-4*(3+4*g2*(-3+4*g2))*g3*g3+8*g1*g1*(1+12*g2*(-1+2*g2)+4*g3-8*g3*g3)+g1*(-3-96*g2*g2+24*g3*(-1+2*g3)+g2*(44+32*(1-2*g3)*g3))));
                    constexpr double w2 = (-1-20*(1+12*g1*(2*g1-1))*w1+20*g3*(2*g3-1)*(4*w1-1))/(20*(1+12*g2*(2*g2-1)+4*g3-8*g3*g3));
                    static const octave_idx_type jk6[][2] = {{1 - 1, 2 - 1}, {1 - 1, 3 - 1}, {1 - 1, 4 - 1}, {2 - 1, 3 - 1}, {2 - 1, 4 - 1}, {3 - 1, 4 - 1}};

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

               return oIntegMass;

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
                    static const double g2[][4] = {{0.5,0.5,0,0},
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

               return oIntegMassDiag;

          default:
               throw std::runtime_error("invalid integration rule");
          }
     }

protected:
     virtual double Jacobian(const ColumnVector& rv, Matrix& J) const {
          FEM_ASSERT(J.rows() == 4);
          FEM_ASSERT(J.columns() == 4);
          FEM_ASSERT(rv.numel() == 4);

          const double Zeta1 = rv.xelem(0);
          const double Zeta2 = rv.xelem(1);
          const double Zeta3 = rv.xelem(2);
          const double Zeta4 = rv.xelem(3);

          FEM_ASSERT(Zeta1 + Zeta2 + Zeta3 + Zeta4 == 1);

          J.xelem(0,0) = 1;
          J.xelem(1,0) = 4*X.xelem(0,7)*Zeta4+4*X.xelem(0,6)*Zeta3+4*X.xelem(0,4)*Zeta2+X.xelem(0,0)*(4*Zeta1-1);
          J.xelem(2,0) = 4*X.xelem(1,7)*Zeta4+4*X.xelem(1,6)*Zeta3+4*X.xelem(1,4)*Zeta2+X.xelem(1,0)*(4*Zeta1-1);
          J.xelem(3,0) = 4*X.xelem(2,7)*Zeta4+4*X.xelem(2,6)*Zeta3+4*X.xelem(2,4)*Zeta2+X.xelem(2,0)*(4*Zeta1-1);
          J.xelem(0,1) = 1;
          J.xelem(1,1) = 4*X.xelem(0,8)*Zeta4+4*X.xelem(0,5)*Zeta3+X.xelem(0,1)*(4*Zeta2-1)+4*X.xelem(0,4)*Zeta1;
          J.xelem(2,1) = 4*X.xelem(1,8)*Zeta4+4*X.xelem(1,5)*Zeta3+X.xelem(1,1)*(4*Zeta2-1)+4*X.xelem(1,4)*Zeta1;
          J.xelem(3,1) = 4*X.xelem(2,8)*Zeta4+4*X.xelem(2,5)*Zeta3+X.xelem(2,1)*(4*Zeta2-1)+4*X.xelem(2,4)*Zeta1;
          J.xelem(0,2) = 1;
          J.xelem(1,2) = 4*X.xelem(0,9)*Zeta4+X.xelem(0,2)*(4*Zeta3-1)+4*X.xelem(0,5)*Zeta2+4*X.xelem(0,6)*Zeta1;
          J.xelem(2,2) = 4*X.xelem(1,9)*Zeta4+X.xelem(1,2)*(4*Zeta3-1)+4*X.xelem(1,5)*Zeta2+4*X.xelem(1,6)*Zeta1;
          J.xelem(3,2) = 4*X.xelem(2,9)*Zeta4+X.xelem(2,2)*(4*Zeta3-1)+4*X.xelem(2,5)*Zeta2+4*X.xelem(2,6)*Zeta1;
          J.xelem(0,3) = 1;
          J.xelem(1,3) = X.xelem(0,3)*(4*Zeta4-1)+4*X.xelem(0,9)*Zeta3+4*X.xelem(0,8)*Zeta2+4*X.xelem(0,7)*Zeta1;
          J.xelem(2,3) = X.xelem(1,3)*(4*Zeta4-1)+4*X.xelem(1,9)*Zeta3+4*X.xelem(1,8)*Zeta2+4*X.xelem(1,7)*Zeta1;
          J.xelem(3,3) = X.xelem(2,3)*(4*Zeta4-1)+4*X.xelem(2,9)*Zeta3+4*X.xelem(2,8)*Zeta2+4*X.xelem(2,7)*Zeta1;

          return J.determinant() * gamma;
     }

     virtual void DispInterpMatrix(const ColumnVector& rv, Matrix& H) const {
          FEM_ASSERT(H.rows() == 3);
          FEM_ASSERT(H.columns() == 30);
          FEM_ASSERT(rv.numel() == 4);

          const double Zeta1 = rv.xelem(0);
          const double Zeta2 = rv.xelem(1);
          const double Zeta3 = rv.xelem(2);
          const double Zeta4 = rv.xelem(3);

          H.xelem(0,0) = Zeta1*(2*Zeta1-1);
          H.xelem(1,0) = 0;
          H.xelem(2,0) = 0;
          H.xelem(0,1) = 0;
          H.xelem(1,1) = Zeta1*(2*Zeta1-1);
          H.xelem(2,1) = 0;
          H.xelem(0,2) = 0;
          H.xelem(1,2) = 0;
          H.xelem(2,2) = Zeta1*(2*Zeta1-1);
          H.xelem(0,3) = Zeta2*(2*Zeta2-1);
          H.xelem(1,3) = 0;
          H.xelem(2,3) = 0;
          H.xelem(0,4) = 0;
          H.xelem(1,4) = Zeta2*(2*Zeta2-1);
          H.xelem(2,4) = 0;
          H.xelem(0,5) = 0;
          H.xelem(1,5) = 0;
          H.xelem(2,5) = Zeta2*(2*Zeta2-1);
          H.xelem(0,6) = Zeta3*(2*Zeta3-1);
          H.xelem(1,6) = 0;
          H.xelem(2,6) = 0;
          H.xelem(0,7) = 0;
          H.xelem(1,7) = Zeta3*(2*Zeta3-1);
          H.xelem(2,7) = 0;
          H.xelem(0,8) = 0;
          H.xelem(1,8) = 0;
          H.xelem(2,8) = Zeta3*(2*Zeta3-1);
          H.xelem(0,9) = Zeta4*(2*Zeta4-1);
          H.xelem(1,9) = 0;
          H.xelem(2,9) = 0;
          H.xelem(0,10) = 0;
          H.xelem(1,10) = Zeta4*(2*Zeta4-1);
          H.xelem(2,10) = 0;
          H.xelem(0,11) = 0;
          H.xelem(1,11) = 0;
          H.xelem(2,11) = Zeta4*(2*Zeta4-1);
          H.xelem(0,12) = 4*Zeta1*Zeta2;
          H.xelem(1,12) = 0;
          H.xelem(2,12) = 0;
          H.xelem(0,13) = 0;
          H.xelem(1,13) = 4*Zeta1*Zeta2;
          H.xelem(2,13) = 0;
          H.xelem(0,14) = 0;
          H.xelem(1,14) = 0;
          H.xelem(2,14) = 4*Zeta1*Zeta2;
          H.xelem(0,15) = 4*Zeta2*Zeta3;
          H.xelem(1,15) = 0;
          H.xelem(2,15) = 0;
          H.xelem(0,16) = 0;
          H.xelem(1,16) = 4*Zeta2*Zeta3;
          H.xelem(2,16) = 0;
          H.xelem(0,17) = 0;
          H.xelem(1,17) = 0;
          H.xelem(2,17) = 4*Zeta2*Zeta3;
          H.xelem(0,18) = 4*Zeta1*Zeta3;
          H.xelem(1,18) = 0;
          H.xelem(2,18) = 0;
          H.xelem(0,19) = 0;
          H.xelem(1,19) = 4*Zeta1*Zeta3;
          H.xelem(2,19) = 0;
          H.xelem(0,20) = 0;
          H.xelem(1,20) = 0;
          H.xelem(2,20) = 4*Zeta1*Zeta3;
          H.xelem(0,21) = 4*Zeta1*Zeta4;
          H.xelem(1,21) = 0;
          H.xelem(2,21) = 0;
          H.xelem(0,22) = 0;
          H.xelem(1,22) = 4*Zeta1*Zeta4;
          H.xelem(2,22) = 0;
          H.xelem(0,23) = 0;
          H.xelem(1,23) = 0;
          H.xelem(2,23) = 4*Zeta1*Zeta4;
          H.xelem(0,24) = 4*Zeta2*Zeta4;
          H.xelem(1,24) = 0;
          H.xelem(2,24) = 0;
          H.xelem(0,25) = 0;
          H.xelem(1,25) = 4*Zeta2*Zeta4;
          H.xelem(2,25) = 0;
          H.xelem(0,26) = 0;
          H.xelem(1,26) = 0;
          H.xelem(2,26) = 4*Zeta2*Zeta4;
          H.xelem(0,27) = 4*Zeta3*Zeta4;
          H.xelem(1,27) = 0;
          H.xelem(2,27) = 0;
          H.xelem(0,28) = 0;
          H.xelem(1,28) = 4*Zeta3*Zeta4;
          H.xelem(2,28) = 0;
          H.xelem(0,29) = 0;
          H.xelem(1,29) = 0;
          H.xelem(2,29) = 4*Zeta3*Zeta4;
     }

     virtual void TemperatureGradientMatrix(const ColumnVector& rv, const Matrix& invJ, Matrix& Bt) const {
          FEM_ASSERT(Bt.rows() == 3);
          FEM_ASSERT(Bt.columns() == 10);
          FEM_ASSERT(rv.numel() == 4);

          const double Zeta1 = rv.xelem(0);
          const double Zeta2 = rv.xelem(1);
          const double Zeta3 = rv.xelem(2);
          const double Zeta4 = rv.xelem(3);
          
          Bt.xelem(0,0) = invJ.xelem(0,1)*(4*Zeta1-1);
          Bt.xelem(1,0) = invJ.xelem(0,2)*(4*Zeta1-1);
          Bt.xelem(2,0) = invJ.xelem(0,3)*(4*Zeta1-1);
          Bt.xelem(0,1) = invJ.xelem(1,1)*(4*Zeta2-1);
          Bt.xelem(1,1) = invJ.xelem(1,2)*(4*Zeta2-1);
          Bt.xelem(2,1) = invJ.xelem(1,3)*(4*Zeta2-1);
          Bt.xelem(0,2) = invJ.xelem(2,1)*(4*Zeta3-1);
          Bt.xelem(1,2) = invJ.xelem(2,2)*(4*Zeta3-1);
          Bt.xelem(2,2) = invJ.xelem(2,3)*(4*Zeta3-1);
          Bt.xelem(0,3) = invJ.xelem(3,1)*(4*Zeta4-1);
          Bt.xelem(1,3) = invJ.xelem(3,2)*(4*Zeta4-1);
          Bt.xelem(2,3) = invJ.xelem(3,3)*(4*Zeta4-1);
          Bt.xelem(0,4) = 4*invJ.xelem(0,1)*Zeta2+4*invJ.xelem(1,1)*Zeta1;
          Bt.xelem(1,4) = 4*invJ.xelem(0,2)*Zeta2+4*invJ.xelem(1,2)*Zeta1;
          Bt.xelem(2,4) = 4*invJ.xelem(0,3)*Zeta2+4*invJ.xelem(1,3)*Zeta1;
          Bt.xelem(0,5) = 4*invJ.xelem(1,1)*Zeta3+4*invJ.xelem(2,1)*Zeta2;
          Bt.xelem(1,5) = 4*invJ.xelem(1,2)*Zeta3+4*invJ.xelem(2,2)*Zeta2;
          Bt.xelem(2,5) = 4*invJ.xelem(1,3)*Zeta3+4*invJ.xelem(2,3)*Zeta2;
          Bt.xelem(0,6) = 4*invJ.xelem(0,1)*Zeta3+4*invJ.xelem(2,1)*Zeta1;
          Bt.xelem(1,6) = 4*invJ.xelem(0,2)*Zeta3+4*invJ.xelem(2,2)*Zeta1;
          Bt.xelem(2,6) = 4*invJ.xelem(0,3)*Zeta3+4*invJ.xelem(2,3)*Zeta1;
          Bt.xelem(0,7) = 4*invJ.xelem(0,1)*Zeta4+4*invJ.xelem(3,1)*Zeta1;
          Bt.xelem(1,7) = 4*invJ.xelem(0,2)*Zeta4+4*invJ.xelem(3,2)*Zeta1;
          Bt.xelem(2,7) = 4*invJ.xelem(0,3)*Zeta4+4*invJ.xelem(3,3)*Zeta1;
          Bt.xelem(0,8) = 4*invJ.xelem(1,1)*Zeta4+4*invJ.xelem(3,1)*Zeta2;
          Bt.xelem(1,8) = 4*invJ.xelem(1,2)*Zeta4+4*invJ.xelem(3,2)*Zeta2;
          Bt.xelem(2,8) = 4*invJ.xelem(1,3)*Zeta4+4*invJ.xelem(3,3)*Zeta2;
          Bt.xelem(0,9) = 4*invJ.xelem(2,1)*Zeta4+4*invJ.xelem(3,1)*Zeta3;
          Bt.xelem(1,9) = 4*invJ.xelem(2,2)*Zeta4+4*invJ.xelem(3,2)*Zeta3;
          Bt.xelem(2,9) = 4*invJ.xelem(2,3)*Zeta4+4*invJ.xelem(3,3)*Zeta3;
     }
     
     virtual void StrainMatrix(const ColumnVector& rv, const Matrix& invJ, Matrix& B) const {
          FEM_ASSERT(invJ.rows() == 4);
          FEM_ASSERT(invJ.columns() == 4);
          FEM_ASSERT(B.rows() == 6);
          FEM_ASSERT(B.columns() == 30);
          FEM_ASSERT(rv.numel() == 4);

          const double Zeta1 = rv.xelem(0);
          const double Zeta2 = rv.xelem(1);
          const double Zeta3 = rv.xelem(2);
          const double Zeta4 = rv.xelem(3);

          B.xelem(0,0) = invJ.xelem(0,1)*(4*Zeta1-1);
          B.xelem(1,0) = 0;
          B.xelem(2,0) = 0;
          B.xelem(3,0) = invJ.xelem(0,2)*(4*Zeta1-1);
          B.xelem(4,0) = 0;
          B.xelem(5,0) = invJ.xelem(0,3)*(4*Zeta1-1);
          B.xelem(0,1) = 0;
          B.xelem(1,1) = invJ.xelem(0,2)*(4*Zeta1-1);
          B.xelem(2,1) = 0;
          B.xelem(3,1) = invJ.xelem(0,1)*(4*Zeta1-1);
          B.xelem(4,1) = invJ.xelem(0,3)*(4*Zeta1-1);
          B.xelem(5,1) = 0;
          B.xelem(0,2) = 0;
          B.xelem(1,2) = 0;
          B.xelem(2,2) = invJ.xelem(0,3)*(4*Zeta1-1);
          B.xelem(3,2) = 0;
          B.xelem(4,2) = invJ.xelem(0,2)*(4*Zeta1-1);
          B.xelem(5,2) = invJ.xelem(0,1)*(4*Zeta1-1);
          B.xelem(0,3) = invJ.xelem(1,1)*(4*Zeta2-1);
          B.xelem(1,3) = 0;
          B.xelem(2,3) = 0;
          B.xelem(3,3) = invJ.xelem(1,2)*(4*Zeta2-1);
          B.xelem(4,3) = 0;
          B.xelem(5,3) = invJ.xelem(1,3)*(4*Zeta2-1);
          B.xelem(0,4) = 0;
          B.xelem(1,4) = invJ.xelem(1,2)*(4*Zeta2-1);
          B.xelem(2,4) = 0;
          B.xelem(3,4) = invJ.xelem(1,1)*(4*Zeta2-1);
          B.xelem(4,4) = invJ.xelem(1,3)*(4*Zeta2-1);
          B.xelem(5,4) = 0;
          B.xelem(0,5) = 0;
          B.xelem(1,5) = 0;
          B.xelem(2,5) = invJ.xelem(1,3)*(4*Zeta2-1);
          B.xelem(3,5) = 0;
          B.xelem(4,5) = invJ.xelem(1,2)*(4*Zeta2-1);
          B.xelem(5,5) = invJ.xelem(1,1)*(4*Zeta2-1);
          B.xelem(0,6) = invJ.xelem(2,1)*(4*Zeta3-1);
          B.xelem(1,6) = 0;
          B.xelem(2,6) = 0;
          B.xelem(3,6) = invJ.xelem(2,2)*(4*Zeta3-1);
          B.xelem(4,6) = 0;
          B.xelem(5,6) = invJ.xelem(2,3)*(4*Zeta3-1);
          B.xelem(0,7) = 0;
          B.xelem(1,7) = invJ.xelem(2,2)*(4*Zeta3-1);
          B.xelem(2,7) = 0;
          B.xelem(3,7) = invJ.xelem(2,1)*(4*Zeta3-1);
          B.xelem(4,7) = invJ.xelem(2,3)*(4*Zeta3-1);
          B.xelem(5,7) = 0;
          B.xelem(0,8) = 0;
          B.xelem(1,8) = 0;
          B.xelem(2,8) = invJ.xelem(2,3)*(4*Zeta3-1);
          B.xelem(3,8) = 0;
          B.xelem(4,8) = invJ.xelem(2,2)*(4*Zeta3-1);
          B.xelem(5,8) = invJ.xelem(2,1)*(4*Zeta3-1);
          B.xelem(0,9) = invJ.xelem(3,1)*(4*Zeta4-1);
          B.xelem(1,9) = 0;
          B.xelem(2,9) = 0;
          B.xelem(3,9) = invJ.xelem(3,2)*(4*Zeta4-1);
          B.xelem(4,9) = 0;
          B.xelem(5,9) = invJ.xelem(3,3)*(4*Zeta4-1);
          B.xelem(0,10) = 0;
          B.xelem(1,10) = invJ.xelem(3,2)*(4*Zeta4-1);
          B.xelem(2,10) = 0;
          B.xelem(3,10) = invJ.xelem(3,1)*(4*Zeta4-1);
          B.xelem(4,10) = invJ.xelem(3,3)*(4*Zeta4-1);
          B.xelem(5,10) = 0;
          B.xelem(0,11) = 0;
          B.xelem(1,11) = 0;
          B.xelem(2,11) = invJ.xelem(3,3)*(4*Zeta4-1);
          B.xelem(3,11) = 0;
          B.xelem(4,11) = invJ.xelem(3,2)*(4*Zeta4-1);
          B.xelem(5,11) = invJ.xelem(3,1)*(4*Zeta4-1);
          B.xelem(0,12) = 4*invJ.xelem(0,1)*Zeta2+4*invJ.xelem(1,1)*Zeta1;
          B.xelem(1,12) = 0;
          B.xelem(2,12) = 0;
          B.xelem(3,12) = 4*invJ.xelem(0,2)*Zeta2+4*invJ.xelem(1,2)*Zeta1;
          B.xelem(4,12) = 0;
          B.xelem(5,12) = 4*invJ.xelem(0,3)*Zeta2+4*invJ.xelem(1,3)*Zeta1;
          B.xelem(0,13) = 0;
          B.xelem(1,13) = 4*invJ.xelem(0,2)*Zeta2+4*invJ.xelem(1,2)*Zeta1;
          B.xelem(2,13) = 0;
          B.xelem(3,13) = 4*invJ.xelem(0,1)*Zeta2+4*invJ.xelem(1,1)*Zeta1;
          B.xelem(4,13) = 4*invJ.xelem(0,3)*Zeta2+4*invJ.xelem(1,3)*Zeta1;
          B.xelem(5,13) = 0;
          B.xelem(0,14) = 0;
          B.xelem(1,14) = 0;
          B.xelem(2,14) = 4*invJ.xelem(0,3)*Zeta2+4*invJ.xelem(1,3)*Zeta1;
          B.xelem(3,14) = 0;
          B.xelem(4,14) = 4*invJ.xelem(0,2)*Zeta2+4*invJ.xelem(1,2)*Zeta1;
          B.xelem(5,14) = 4*invJ.xelem(0,1)*Zeta2+4*invJ.xelem(1,1)*Zeta1;
          B.xelem(0,15) = 4*invJ.xelem(1,1)*Zeta3+4*invJ.xelem(2,1)*Zeta2;
          B.xelem(1,15) = 0;
          B.xelem(2,15) = 0;
          B.xelem(3,15) = 4*invJ.xelem(1,2)*Zeta3+4*invJ.xelem(2,2)*Zeta2;
          B.xelem(4,15) = 0;
          B.xelem(5,15) = 4*invJ.xelem(1,3)*Zeta3+4*invJ.xelem(2,3)*Zeta2;
          B.xelem(0,16) = 0;
          B.xelem(1,16) = 4*invJ.xelem(1,2)*Zeta3+4*invJ.xelem(2,2)*Zeta2;
          B.xelem(2,16) = 0;
          B.xelem(3,16) = 4*invJ.xelem(1,1)*Zeta3+4*invJ.xelem(2,1)*Zeta2;
          B.xelem(4,16) = 4*invJ.xelem(1,3)*Zeta3+4*invJ.xelem(2,3)*Zeta2;
          B.xelem(5,16) = 0;
          B.xelem(0,17) = 0;
          B.xelem(1,17) = 0;
          B.xelem(2,17) = 4*invJ.xelem(1,3)*Zeta3+4*invJ.xelem(2,3)*Zeta2;
          B.xelem(3,17) = 0;
          B.xelem(4,17) = 4*invJ.xelem(1,2)*Zeta3+4*invJ.xelem(2,2)*Zeta2;
          B.xelem(5,17) = 4*invJ.xelem(1,1)*Zeta3+4*invJ.xelem(2,1)*Zeta2;
          B.xelem(0,18) = 4*invJ.xelem(0,1)*Zeta3+4*invJ.xelem(2,1)*Zeta1;
          B.xelem(1,18) = 0;
          B.xelem(2,18) = 0;
          B.xelem(3,18) = 4*invJ.xelem(0,2)*Zeta3+4*invJ.xelem(2,2)*Zeta1;
          B.xelem(4,18) = 0;
          B.xelem(5,18) = 4*invJ.xelem(0,3)*Zeta3+4*invJ.xelem(2,3)*Zeta1;
          B.xelem(0,19) = 0;
          B.xelem(1,19) = 4*invJ.xelem(0,2)*Zeta3+4*invJ.xelem(2,2)*Zeta1;
          B.xelem(2,19) = 0;
          B.xelem(3,19) = 4*invJ.xelem(0,1)*Zeta3+4*invJ.xelem(2,1)*Zeta1;
          B.xelem(4,19) = 4*invJ.xelem(0,3)*Zeta3+4*invJ.xelem(2,3)*Zeta1;
          B.xelem(5,19) = 0;
          B.xelem(0,20) = 0;
          B.xelem(1,20) = 0;
          B.xelem(2,20) = 4*invJ.xelem(0,3)*Zeta3+4*invJ.xelem(2,3)*Zeta1;
          B.xelem(3,20) = 0;
          B.xelem(4,20) = 4*invJ.xelem(0,2)*Zeta3+4*invJ.xelem(2,2)*Zeta1;
          B.xelem(5,20) = 4*invJ.xelem(0,1)*Zeta3+4*invJ.xelem(2,1)*Zeta1;
          B.xelem(0,21) = 4*invJ.xelem(0,1)*Zeta4+4*invJ.xelem(3,1)*Zeta1;
          B.xelem(1,21) = 0;
          B.xelem(2,21) = 0;
          B.xelem(3,21) = 4*invJ.xelem(0,2)*Zeta4+4*invJ.xelem(3,2)*Zeta1;
          B.xelem(4,21) = 0;
          B.xelem(5,21) = 4*invJ.xelem(0,3)*Zeta4+4*invJ.xelem(3,3)*Zeta1;
          B.xelem(0,22) = 0;
          B.xelem(1,22) = 4*invJ.xelem(0,2)*Zeta4+4*invJ.xelem(3,2)*Zeta1;
          B.xelem(2,22) = 0;
          B.xelem(3,22) = 4*invJ.xelem(0,1)*Zeta4+4*invJ.xelem(3,1)*Zeta1;
          B.xelem(4,22) = 4*invJ.xelem(0,3)*Zeta4+4*invJ.xelem(3,3)*Zeta1;
          B.xelem(5,22) = 0;
          B.xelem(0,23) = 0;
          B.xelem(1,23) = 0;
          B.xelem(2,23) = 4*invJ.xelem(0,3)*Zeta4+4*invJ.xelem(3,3)*Zeta1;
          B.xelem(3,23) = 0;
          B.xelem(4,23) = 4*invJ.xelem(0,2)*Zeta4+4*invJ.xelem(3,2)*Zeta1;
          B.xelem(5,23) = 4*invJ.xelem(0,1)*Zeta4+4*invJ.xelem(3,1)*Zeta1;
          B.xelem(0,24) = 4*invJ.xelem(1,1)*Zeta4+4*invJ.xelem(3,1)*Zeta2;
          B.xelem(1,24) = 0;
          B.xelem(2,24) = 0;
          B.xelem(3,24) = 4*invJ.xelem(1,2)*Zeta4+4*invJ.xelem(3,2)*Zeta2;
          B.xelem(4,24) = 0;
          B.xelem(5,24) = 4*invJ.xelem(1,3)*Zeta4+4*invJ.xelem(3,3)*Zeta2;
          B.xelem(0,25) = 0;
          B.xelem(1,25) = 4*invJ.xelem(1,2)*Zeta4+4*invJ.xelem(3,2)*Zeta2;
          B.xelem(2,25) = 0;
          B.xelem(3,25) = 4*invJ.xelem(1,1)*Zeta4+4*invJ.xelem(3,1)*Zeta2;
          B.xelem(4,25) = 4*invJ.xelem(1,3)*Zeta4+4*invJ.xelem(3,3)*Zeta2;
          B.xelem(5,25) = 0;
          B.xelem(0,26) = 0;
          B.xelem(1,26) = 0;
          B.xelem(2,26) = 4*invJ.xelem(1,3)*Zeta4+4*invJ.xelem(3,3)*Zeta2;
          B.xelem(3,26) = 0;
          B.xelem(4,26) = 4*invJ.xelem(1,2)*Zeta4+4*invJ.xelem(3,2)*Zeta2;
          B.xelem(5,26) = 4*invJ.xelem(1,1)*Zeta4+4*invJ.xelem(3,1)*Zeta2;
          B.xelem(0,27) = 4*invJ.xelem(2,1)*Zeta4+4*invJ.xelem(3,1)*Zeta3;
          B.xelem(1,27) = 0;
          B.xelem(2,27) = 0;
          B.xelem(3,27) = 4*invJ.xelem(2,2)*Zeta4+4*invJ.xelem(3,2)*Zeta3;
          B.xelem(4,27) = 0;
          B.xelem(5,27) = 4*invJ.xelem(2,3)*Zeta4+4*invJ.xelem(3,3)*Zeta3;
          B.xelem(0,28) = 0;
          B.xelem(1,28) = 4*invJ.xelem(2,2)*Zeta4+4*invJ.xelem(3,2)*Zeta3;
          B.xelem(2,28) = 0;
          B.xelem(3,28) = 4*invJ.xelem(2,1)*Zeta4+4*invJ.xelem(3,1)*Zeta3;
          B.xelem(4,28) = 4*invJ.xelem(2,3)*Zeta4+4*invJ.xelem(3,3)*Zeta3;
          B.xelem(5,28) = 0;
          B.xelem(0,29) = 0;
          B.xelem(1,29) = 0;
          B.xelem(2,29) = 4*invJ.xelem(2,3)*Zeta4+4*invJ.xelem(3,3)*Zeta3;
          B.xelem(3,29) = 0;
          B.xelem(4,29) = 4*invJ.xelem(2,2)*Zeta4+4*invJ.xelem(3,2)*Zeta3;
          B.xelem(5,29) = 4*invJ.xelem(2,1)*Zeta4+4*invJ.xelem(3,1)*Zeta3;
     }

     virtual Matrix InterpGaussToNodal(FemMatrixType eMatType, const Matrix& taug) const {
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

          Matrix taun = H.solve(taug);

          taun.resize(nodes.numel(), taun.columns());

          static const struct {
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

     void ScalarInterpMatrix(const ColumnVector& rv, Matrix& Hs, octave_idx_type irow) const {
          FEM_ASSERT(rv.numel() == 4);
          FEM_ASSERT(Hs.columns() == 10);
          FEM_ASSERT(irow >= 0);
          FEM_ASSERT(irow < Hs.rows());

          const double Zeta1 = rv.xelem(0);
          const double Zeta2 = rv.xelem(1);
          const double Zeta3 = rv.xelem(2);
          const double Zeta4 = rv.xelem(3);

          Hs.xelem(irow, 0) = Zeta1*(2*Zeta1-1);
          Hs.xelem(irow, 1) = Zeta2*(2*Zeta2-1);
          Hs.xelem(irow, 2) = Zeta3*(2*Zeta3-1);
          Hs.xelem(irow, 3) = Zeta4*(2*Zeta4-1);
          Hs.xelem(irow, 4) = 4*Zeta1*Zeta2;
          Hs.xelem(irow, 5) = 4*Zeta2*Zeta3;
          Hs.xelem(irow, 6) = 4*Zeta1*Zeta3;
          Hs.xelem(irow, 7) = 4*Zeta1*Zeta4;
          Hs.xelem(irow, 8) = 4*Zeta2*Zeta4;
          Hs.xelem(irow, 9) = 4*Zeta3*Zeta4;
     }
     
private:
     void ScalarInterpMatrixCornerNodes(const ColumnVector& rv, Matrix& Hs, octave_idx_type irow) const {
          FEM_ASSERT(rv.numel() == 4);
          FEM_ASSERT(Hs.columns() == 4);
          FEM_ASSERT(irow >= 0);
          FEM_ASSERT(irow < Hs.rows());

          const double Zeta1 = rv.xelem(0);
          const double Zeta2 = rv.xelem(1);
          const double Zeta3 = rv.xelem(2);
          const double Zeta4 = rv.xelem(3);

          Hs.xelem(irow, 0) = Zeta1;
          Hs.xelem(irow, 1) = Zeta2;
          Hs.xelem(irow, 2) = Zeta3;
          Hs.xelem(irow, 3) = Zeta4;
     }
};

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

     static double EqualityConstr(const ColumnVector& rv) {
          return rv.xelem(0) + rv.xelem(1) + rv.xelem(2) - 1;
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

     static void ScalarInterpMatrix(const ColumnVector& rv, RowVector& HA) {
          FEM_ASSERT(rv.rows() == 3);
          FEM_ASSERT(HA.columns() == 6);

          const double Zeta1 = rv.xelem(0);
          const double Zeta2 = rv.xelem(1);
          const double Zeta3 = rv.xelem(2);

          HA.xelem(0) = Zeta1*(2*Zeta1-1);
          HA.xelem(1) = Zeta2*(2*Zeta2-1);
          HA.xelem(2) = Zeta3*(2*Zeta3-1);
          HA.xelem(3) = 4*Zeta1*Zeta2;
          HA.xelem(4) = 4*Zeta2*Zeta3;
          HA.xelem(5) = 4*Zeta1*Zeta3;
     }

     static void VectorInterpMatrix(const ColumnVector& rv, Matrix& Hf) {
          FEM_ASSERT(rv.rows() == 3);
          FEM_ASSERT(Hf.rows() == 3);
          FEM_ASSERT(Hf.columns() == 18);

          const double Zeta1 = rv.xelem(0);
          const double Zeta2 = rv.xelem(1);
          const double Zeta3 = rv.xelem(2);

          Hf.xelem(0,0) = Zeta1*(2*Zeta1-1);
          Hf.xelem(1,0) = 0;
          Hf.xelem(2,0) = 0;
          Hf.xelem(0,1) = 0;
          Hf.xelem(1,1) = Zeta1*(2*Zeta1-1);
          Hf.xelem(2,1) = 0;
          Hf.xelem(0,2) = 0;
          Hf.xelem(1,2) = 0;
          Hf.xelem(2,2) = Zeta1*(2*Zeta1-1);
          Hf.xelem(0,3) = Zeta2*(2*Zeta2-1);
          Hf.xelem(1,3) = 0;
          Hf.xelem(2,3) = 0;
          Hf.xelem(0,4) = 0;
          Hf.xelem(1,4) = Zeta2*(2*Zeta2-1);
          Hf.xelem(2,4) = 0;
          Hf.xelem(0,5) = 0;
          Hf.xelem(1,5) = 0;
          Hf.xelem(2,5) = Zeta2*(2*Zeta2-1);
          Hf.xelem(0,6) = Zeta3*(2*Zeta3-1);
          Hf.xelem(1,6) = 0;
          Hf.xelem(2,6) = 0;
          Hf.xelem(0,7) = 0;
          Hf.xelem(1,7) = Zeta3*(2*Zeta3-1);
          Hf.xelem(2,7) = 0;
          Hf.xelem(0,8) = 0;
          Hf.xelem(1,8) = 0;
          Hf.xelem(2,8) = Zeta3*(2*Zeta3-1);
          Hf.xelem(0,9) = 4*Zeta1*Zeta2;
          Hf.xelem(1,9) = 0;
          Hf.xelem(2,9) = 0;
          Hf.xelem(0,10) = 0;
          Hf.xelem(1,10) = 4*Zeta1*Zeta2;
          Hf.xelem(2,10) = 0;
          Hf.xelem(0,11) = 0;
          Hf.xelem(1,11) = 0;
          Hf.xelem(2,11) = 4*Zeta1*Zeta2;
          Hf.xelem(0,12) = 4*Zeta2*Zeta3;
          Hf.xelem(1,12) = 0;
          Hf.xelem(2,12) = 0;
          Hf.xelem(0,13) = 0;
          Hf.xelem(1,13) = 4*Zeta2*Zeta3;
          Hf.xelem(2,13) = 0;
          Hf.xelem(0,14) = 0;
          Hf.xelem(1,14) = 0;
          Hf.xelem(2,14) = 4*Zeta2*Zeta3;
          Hf.xelem(0,15) = 4*Zeta1*Zeta3;
          Hf.xelem(1,15) = 0;
          Hf.xelem(2,15) = 0;
          Hf.xelem(0,16) = 0;
          Hf.xelem(1,16) = 4*Zeta1*Zeta3;
          Hf.xelem(2,16) = 0;
          Hf.xelem(0,17) = 0;
          Hf.xelem(1,17) = 0;
          Hf.xelem(2,17) = 4*Zeta1*Zeta3;
     }

     static void VectorInterpMatrixDerR(const ColumnVector& rv, Matrix& dHf_dr) {
          FEM_ASSERT(rv.rows() == 3);
          FEM_ASSERT(dHf_dr.rows() == 3);
          FEM_ASSERT(dHf_dr.columns() == 18);

          const double Zeta1 = rv.xelem(0);
          const double Zeta2 = rv.xelem(1);
          const double Zeta3 = rv.xelem(2);

          dHf_dr.xelem(0,0) = 4*Zeta1-1;
          dHf_dr.xelem(1,0) = 0;
          dHf_dr.xelem(2,0) = 0;
          dHf_dr.xelem(0,1) = 0;
          dHf_dr.xelem(1,1) = 4*Zeta1-1;
          dHf_dr.xelem(2,1) = 0;
          dHf_dr.xelem(0,2) = 0;
          dHf_dr.xelem(1,2) = 0;
          dHf_dr.xelem(2,2) = 4*Zeta1-1;
          dHf_dr.xelem(0,3) = 0;
          dHf_dr.xelem(1,3) = 0;
          dHf_dr.xelem(2,3) = 0;
          dHf_dr.xelem(0,4) = 0;
          dHf_dr.xelem(1,4) = 0;
          dHf_dr.xelem(2,4) = 0;
          dHf_dr.xelem(0,5) = 0;
          dHf_dr.xelem(1,5) = 0;
          dHf_dr.xelem(2,5) = 0;
          dHf_dr.xelem(0,6) = 1-4*Zeta3;
          dHf_dr.xelem(1,6) = 0;
          dHf_dr.xelem(2,6) = 0;
          dHf_dr.xelem(0,7) = 0;
          dHf_dr.xelem(1,7) = 1-4*Zeta3;
          dHf_dr.xelem(2,7) = 0;
          dHf_dr.xelem(0,8) = 0;
          dHf_dr.xelem(1,8) = 0;
          dHf_dr.xelem(2,8) = 1-4*Zeta3;
          dHf_dr.xelem(0,9) = 4*Zeta2;
          dHf_dr.xelem(1,9) = 0;
          dHf_dr.xelem(2,9) = 0;
          dHf_dr.xelem(0,10) = 0;
          dHf_dr.xelem(1,10) = 4*Zeta2;
          dHf_dr.xelem(2,10) = 0;
          dHf_dr.xelem(0,11) = 0;
          dHf_dr.xelem(1,11) = 0;
          dHf_dr.xelem(2,11) = 4*Zeta2;
          dHf_dr.xelem(0,12) = -4*Zeta2;
          dHf_dr.xelem(1,12) = 0;
          dHf_dr.xelem(2,12) = 0;
          dHf_dr.xelem(0,13) = 0;
          dHf_dr.xelem(1,13) = -4*Zeta2;
          dHf_dr.xelem(2,13) = 0;
          dHf_dr.xelem(0,14) = 0;
          dHf_dr.xelem(1,14) = 0;
          dHf_dr.xelem(2,14) = -4*Zeta2;
          dHf_dr.xelem(0,15) = 4*Zeta3-4*Zeta1;
          dHf_dr.xelem(1,15) = 0;
          dHf_dr.xelem(2,15) = 0;
          dHf_dr.xelem(0,16) = 0;
          dHf_dr.xelem(1,16) = 4*Zeta3-4*Zeta1;
          dHf_dr.xelem(2,16) = 0;
          dHf_dr.xelem(0,17) = 0;
          dHf_dr.xelem(1,17) = 0;
          dHf_dr.xelem(2,17) = 4*Zeta3-4*Zeta1;
     }

     static void VectorInterpMatrixDerS(const ColumnVector& rv, Matrix& dHf_ds) {
          FEM_ASSERT(rv.rows() == 3);
          FEM_ASSERT(dHf_ds.rows() == 3);
          FEM_ASSERT(dHf_ds.columns() == 18);

          const double Zeta1 = rv.xelem(0);
          const double Zeta2 = rv.xelem(1);
          const double Zeta3 = rv.xelem(2);

          dHf_ds.xelem(0,0) = 0;
          dHf_ds.xelem(1,0) = 0;
          dHf_ds.xelem(2,0) = 0;
          dHf_ds.xelem(0,1) = 0;
          dHf_ds.xelem(1,1) = 0;
          dHf_ds.xelem(2,1) = 0;
          dHf_ds.xelem(0,2) = 0;
          dHf_ds.xelem(1,2) = 0;
          dHf_ds.xelem(2,2) = 0;
          dHf_ds.xelem(0,3) = 4*Zeta2-1;
          dHf_ds.xelem(1,3) = 0;
          dHf_ds.xelem(2,3) = 0;
          dHf_ds.xelem(0,4) = 0;
          dHf_ds.xelem(1,4) = 4*Zeta2-1;
          dHf_ds.xelem(2,4) = 0;
          dHf_ds.xelem(0,5) = 0;
          dHf_ds.xelem(1,5) = 0;
          dHf_ds.xelem(2,5) = 4*Zeta2-1;
          dHf_ds.xelem(0,6) = 1-4*Zeta3;
          dHf_ds.xelem(1,6) = 0;
          dHf_ds.xelem(2,6) = 0;
          dHf_ds.xelem(0,7) = 0;
          dHf_ds.xelem(1,7) = 1-4*Zeta3;
          dHf_ds.xelem(2,7) = 0;
          dHf_ds.xelem(0,8) = 0;
          dHf_ds.xelem(1,8) = 0;
          dHf_ds.xelem(2,8) = 1-4*Zeta3;
          dHf_ds.xelem(0,9) = 4*Zeta1;
          dHf_ds.xelem(1,9) = 0;
          dHf_ds.xelem(2,9) = 0;
          dHf_ds.xelem(0,10) = 0;
          dHf_ds.xelem(1,10) = 4*Zeta1;
          dHf_ds.xelem(2,10) = 0;
          dHf_ds.xelem(0,11) = 0;
          dHf_ds.xelem(1,11) = 0;
          dHf_ds.xelem(2,11) = 4*Zeta1;
          dHf_ds.xelem(0,12) = 4*Zeta3-4*Zeta2;
          dHf_ds.xelem(1,12) = 0;
          dHf_ds.xelem(2,12) = 0;
          dHf_ds.xelem(0,13) = 0;
          dHf_ds.xelem(1,13) = 4*Zeta3-4*Zeta2;
          dHf_ds.xelem(2,13) = 0;
          dHf_ds.xelem(0,14) = 0;
          dHf_ds.xelem(1,14) = 0;
          dHf_ds.xelem(2,14) = 4*Zeta3-4*Zeta2;
          dHf_ds.xelem(0,15) = -4*Zeta1;
          dHf_ds.xelem(1,15) = 0;
          dHf_ds.xelem(2,15) = 0;
          dHf_ds.xelem(0,16) = 0;
          dHf_ds.xelem(1,16) = -4*Zeta1;
          dHf_ds.xelem(2,16) = 0;
          dHf_ds.xelem(0,17) = 0;
          dHf_ds.xelem(1,17) = 0;
          dHf_ds.xelem(2,17) = -4*Zeta1;
     }

     static const IntegrationRule& GetIntegrationRule(Element::FemMatrixType eMatType) {
          static IntegrationRule oIntegLumped, oIntegConsistent;
          constexpr double tria_area = 0.5; // Factor for triangular area

          switch (eMatType) {
          case Element::VEC_LOAD_LUMPED: {
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
                    static const double g2[][3] = {{0.5, 0.5, 0.0},
                                                   {0.0, 0.5, 0.5},
                                                   {0.5, 0.0, 0.5}};

                    for (octave_idx_type i = 0; i < 3; ++i) {
                         oIntegLumped.SetWeight(i + 3, tria_area * w2); // Factor 0.5 for triangle

                         for (octave_idx_type j = 0; j < 3; ++j) {
                              oIntegLumped.SetPosition(i + 3, j, g2[i][j]);
                         }
                    }
               }

               return oIntegLumped;
          } break;
          case Element::VEC_LOAD_CONSISTENT:
          case Element::VEC_LOAD_THERMAL:
          case Element::MAT_THERMAL_COND: {
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

               return oIntegConsistent;
          } break;
          default:
               throw std::runtime_error("invalid matrix type");
          }
     }
};

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

     static constexpr double EqualityConstr(const ColumnVector&) {
          return 0;
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

     static void ScalarInterpMatrix(const ColumnVector& rv, RowVector& HA) {
          FEM_ASSERT(rv.rows() == 2);

          const double r = rv.xelem(0);
          const double s = rv.xelem(1);

          HA.xelem(0) = ((r+1)*(s+1))/4.0;
          HA.xelem(1) = ((1-r)*(s+1))/4.0;
          HA.xelem(2) = ((1-r)*(1-s))/4.0;
          HA.xelem(3) = ((r+1)*(1-s))/4.0;
     }

     static void VectorInterpMatrix(const ColumnVector& rv, Matrix& Hf) {
          FEM_ASSERT(rv.rows() == 2);
          FEM_ASSERT(Hf.rows() == 3);
          FEM_ASSERT(Hf.columns() == 12);

          const double r = rv.xelem(0);
          const double s = rv.xelem(1);

          Hf.xelem(0,0) = ((r+1)*(s+1))/4.0;
          Hf.xelem(1,0) = 0;
          Hf.xelem(2,0) = 0;
          Hf.xelem(0,1) = 0;
          Hf.xelem(1,1) = ((r+1)*(s+1))/4.0;
          Hf.xelem(2,1) = 0;
          Hf.xelem(0,2) = 0;
          Hf.xelem(1,2) = 0;
          Hf.xelem(2,2) = ((r+1)*(s+1))/4.0;
          Hf.xelem(0,3) = ((1-r)*(s+1))/4.0;
          Hf.xelem(1,3) = 0;
          Hf.xelem(2,3) = 0;
          Hf.xelem(0,4) = 0;
          Hf.xelem(1,4) = ((1-r)*(s+1))/4.0;
          Hf.xelem(2,4) = 0;
          Hf.xelem(0,5) = 0;
          Hf.xelem(1,5) = 0;
          Hf.xelem(2,5) = ((1-r)*(s+1))/4.0;
          Hf.xelem(0,6) = ((1-r)*(1-s))/4.0;
          Hf.xelem(1,6) = 0;
          Hf.xelem(2,6) = 0;
          Hf.xelem(0,7) = 0;
          Hf.xelem(1,7) = ((1-r)*(1-s))/4.0;
          Hf.xelem(2,7) = 0;
          Hf.xelem(0,8) = 0;
          Hf.xelem(1,8) = 0;
          Hf.xelem(2,8) = ((1-r)*(1-s))/4.0;
          Hf.xelem(0,9) = ((r+1)*(1-s))/4.0;
          Hf.xelem(1,9) = 0;
          Hf.xelem(2,9) = 0;
          Hf.xelem(0,10) = 0;
          Hf.xelem(1,10) = ((r+1)*(1-s))/4.0;
          Hf.xelem(2,10) = 0;
          Hf.xelem(0,11) = 0;
          Hf.xelem(1,11) = 0;
          Hf.xelem(2,11) = ((r+1)*(1-s))/4.0;

     }

     static void VectorInterpMatrixDerR(const ColumnVector& rv, Matrix& dHf_dr) {
          FEM_ASSERT(rv.rows() == 2);

          // const double r = rv(0);
          const double s = rv.xelem(1);

          dHf_dr.xelem(0,0) = (s+1)/4.0;
          dHf_dr.xelem(1,0) = 0;
          dHf_dr.xelem(2,0) = 0;
          dHf_dr.xelem(0,1) = 0;
          dHf_dr.xelem(1,1) = (s+1)/4.0;
          dHf_dr.xelem(2,1) = 0;
          dHf_dr.xelem(0,2) = 0;
          dHf_dr.xelem(1,2) = 0;
          dHf_dr.xelem(2,2) = (s+1)/4.0;
          dHf_dr.xelem(0,3) = -(s+1)/4.0;
          dHf_dr.xelem(1,3) = 0;
          dHf_dr.xelem(2,3) = 0;
          dHf_dr.xelem(0,4) = 0;
          dHf_dr.xelem(1,4) = -(s+1)/4.0;
          dHf_dr.xelem(2,4) = 0;
          dHf_dr.xelem(0,5) = 0;
          dHf_dr.xelem(1,5) = 0;
          dHf_dr.xelem(2,5) = -(s+1)/4.0;
          dHf_dr.xelem(0,6) = -(1-s)/4.0;
          dHf_dr.xelem(1,6) = 0;
          dHf_dr.xelem(2,6) = 0;
          dHf_dr.xelem(0,7) = 0;
          dHf_dr.xelem(1,7) = -(1-s)/4.0;
          dHf_dr.xelem(2,7) = 0;
          dHf_dr.xelem(0,8) = 0;
          dHf_dr.xelem(1,8) = 0;
          dHf_dr.xelem(2,8) = -(1-s)/4.0;
          dHf_dr.xelem(0,9) = (1-s)/4.0;
          dHf_dr.xelem(1,9) = 0;
          dHf_dr.xelem(2,9) = 0;
          dHf_dr.xelem(0,10) = 0;
          dHf_dr.xelem(1,10) = (1-s)/4.0;
          dHf_dr.xelem(2,10) = 0;
          dHf_dr.xelem(0,11) = 0;
          dHf_dr.xelem(1,11) = 0;
          dHf_dr.xelem(2,11) = (1-s)/4.0;
     }

     static void VectorInterpMatrixDerS(const ColumnVector& rv, Matrix& dHf_ds) {
          FEM_ASSERT(rv.rows() == 2);

          const double r = rv.xelem(0);
          // const double s = rv(1);

          dHf_ds.xelem(0,0) = (r+1)/4.0;
          dHf_ds.xelem(1,0) = 0;
          dHf_ds.xelem(2,0) = 0;
          dHf_ds.xelem(0,1) = 0;
          dHf_ds.xelem(1,1) = (r+1)/4.0;
          dHf_ds.xelem(2,1) = 0;
          dHf_ds.xelem(0,2) = 0;
          dHf_ds.xelem(1,2) = 0;
          dHf_ds.xelem(2,2) = (r+1)/4.0;
          dHf_ds.xelem(0,3) = (1-r)/4.0;
          dHf_ds.xelem(1,3) = 0;
          dHf_ds.xelem(2,3) = 0;
          dHf_ds.xelem(0,4) = 0;
          dHf_ds.xelem(1,4) = (1-r)/4.0;
          dHf_ds.xelem(2,4) = 0;
          dHf_ds.xelem(0,5) = 0;
          dHf_ds.xelem(1,5) = 0;
          dHf_ds.xelem(2,5) = (1-r)/4.0;
          dHf_ds.xelem(0,6) = -(1-r)/4.0;
          dHf_ds.xelem(1,6) = 0;
          dHf_ds.xelem(2,6) = 0;
          dHf_ds.xelem(0,7) = 0;
          dHf_ds.xelem(1,7) = -(1-r)/4.0;
          dHf_ds.xelem(2,7) = 0;
          dHf_ds.xelem(0,8) = 0;
          dHf_ds.xelem(1,8) = 0;
          dHf_ds.xelem(2,8) = -(1-r)/4.0;
          dHf_ds.xelem(0,9) = -(r+1)/4.0;
          dHf_ds.xelem(1,9) = 0;
          dHf_ds.xelem(2,9) = 0;
          dHf_ds.xelem(0,10) = 0;
          dHf_ds.xelem(1,10) = -(r+1)/4.0;
          dHf_ds.xelem(2,10) = 0;
          dHf_ds.xelem(0,11) = 0;
          dHf_ds.xelem(1,11) = 0;
          dHf_ds.xelem(2,11) = -(r+1)/4.0;
     }

     static const IntegrationRule& GetIntegrationRule(Element::FemMatrixType eMatType) {
          static const octave_idx_type N = 2;
          static const double r[2][N] = {{0.577350269189626, -0.577350269189626}, {1., -1.}};
          static const double alpha[2][N] = {{1., 1.}, {1., 1.}};

          static array<IntegrationRule, 2> rgIntegRule;

          octave_idx_type iIntegRule;

          switch (eMatType) {
          case Element::VEC_LOAD_LUMPED:
               iIntegRule = 1;
               break;
          default:
               iIntegRule = 0;
          }

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

          return rgIntegRule[iIntegRule];
     }
};

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

     static constexpr double EqualityConstr(const ColumnVector&) {
          return 0;
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

     static void ScalarInterpMatrix(const ColumnVector& rv, RowVector& HA) {
          FEM_ASSERT(rv.rows() == 2);

          const double r = rv.xelem(0);
          const double s = rv.xelem(1);
          const double r2 = r * r;
          const double s2 = s * s;

          HA.xelem(0) = ((r+1)*(s+1))/4.0E+0-(((r+1)*(1-s2))/2.0E+0+((1-r2)*(s+1))/2.0E+0)/2.0E+0;
          HA.xelem(1) = ((1-r)*(s+1))/4.0E+0-(((1-r)*(1-s2))/2.0E+0+((1-r2)*(s+1))/2.0E+0)/2.0E+0;
          HA.xelem(2) = ((1-r)*(1-s))/4.0E+0-(((1-r)*(1-s2))/2.0E+0+((1-r2)*(1-s))/2.0E+0)/2.0E+0;
          HA.xelem(3) = ((r+1)*(1-s))/4.0E+0-(((r+1)*(1-s2))/2.0E+0+((1-r2)*(1-s))/2.0E+0)/2.0E+0;
          HA.xelem(4) = ((1-r2)*(s+1))/2.0E+0;
          HA.xelem(5) = ((1-r)*(1-s2))/2.0E+0;
          HA.xelem(6) = ((1-r2)*(1-s))/2.0E+0;
          HA.xelem(7) = ((r+1)*(1-s2))/2.0E+0;
     }

     static void VectorInterpMatrix(const ColumnVector& rv, Matrix& Hf) {
          FEM_ASSERT(rv.rows() == 2);

          const double r = rv.xelem(0);
          const double s = rv.xelem(1);
          const double r2 = r * r;
          const double s2 = s * s;

          Hf.xelem(0,0) = ((r+1)*(s+1))/4.0E+0-(((r+1)*(1-s2))/2.0E+0+((1-r2)*(s+1))/2.0E+0)/2.0E+0;
          Hf.xelem(1,0) = 0;
          Hf.xelem(2,0) = 0;
          Hf.xelem(0,1) = 0;
          Hf.xelem(1,1) = ((r+1)*(s+1))/4.0E+0-(((r+1)*(1-s2))/2.0E+0+((1-r2)*(s+1))/2.0E+0)/2.0E+0;
          Hf.xelem(2,1) = 0;
          Hf.xelem(0,2) = 0;
          Hf.xelem(1,2) = 0;
          Hf.xelem(2,2) = ((r+1)*(s+1))/4.0E+0-(((r+1)*(1-s2))/2.0E+0+((1-r2)*(s+1))/2.0E+0)/2.0E+0;
          Hf.xelem(0,3) = ((1-r)*(s+1))/4.0E+0-(((1-r)*(1-s2))/2.0E+0+((1-r2)*(s+1))/2.0E+0)/2.0E+0;
          Hf.xelem(1,3) = 0;
          Hf.xelem(2,3) = 0;
          Hf.xelem(0,4) = 0;
          Hf.xelem(1,4) = ((1-r)*(s+1))/4.0E+0-(((1-r)*(1-s2))/2.0E+0+((1-r2)*(s+1))/2.0E+0)/2.0E+0;
          Hf.xelem(2,4) = 0;
          Hf.xelem(0,5) = 0;
          Hf.xelem(1,5) = 0;
          Hf.xelem(2,5) = ((1-r)*(s+1))/4.0E+0-(((1-r)*(1-s2))/2.0E+0+((1-r2)*(s+1))/2.0E+0)/2.0E+0;
          Hf.xelem(0,6) = ((1-r)*(1-s))/4.0E+0-(((1-r)*(1-s2))/2.0E+0+((1-r2)*(1-s))/2.0E+0)/2.0E+0;
          Hf.xelem(1,6) = 0;
          Hf.xelem(2,6) = 0;
          Hf.xelem(0,7) = 0;
          Hf.xelem(1,7) = ((1-r)*(1-s))/4.0E+0-(((1-r)*(1-s2))/2.0E+0+((1-r2)*(1-s))/2.0E+0)/2.0E+0;
          Hf.xelem(2,7) = 0;
          Hf.xelem(0,8) = 0;
          Hf.xelem(1,8) = 0;
          Hf.xelem(2,8) = ((1-r)*(1-s))/4.0E+0-(((1-r)*(1-s2))/2.0E+0+((1-r2)*(1-s))/2.0E+0)/2.0E+0;
          Hf.xelem(0,9) = ((r+1)*(1-s))/4.0E+0-(((r+1)*(1-s2))/2.0E+0+((1-r2)*(1-s))/2.0E+0)/2.0E+0;
          Hf.xelem(1,9) = 0;
          Hf.xelem(2,9) = 0;
          Hf.xelem(0,10) = 0;
          Hf.xelem(1,10) = ((r+1)*(1-s))/4.0E+0-(((r+1)*(1-s2))/2.0E+0+((1-r2)*(1-s))/2.0E+0)/2.0E+0;
          Hf.xelem(2,10) = 0;
          Hf.xelem(0,11) = 0;
          Hf.xelem(1,11) = 0;
          Hf.xelem(2,11) = ((r+1)*(1-s))/4.0E+0-(((r+1)*(1-s2))/2.0E+0+((1-r2)*(1-s))/2.0E+0)/2.0E+0;
          Hf.xelem(0,12) = ((1-r2)*(s+1))/2.0E+0;
          Hf.xelem(1,12) = 0;
          Hf.xelem(2,12) = 0;
          Hf.xelem(0,13) = 0;
          Hf.xelem(1,13) = ((1-r2)*(s+1))/2.0E+0;
          Hf.xelem(2,13) = 0;
          Hf.xelem(0,14) = 0;
          Hf.xelem(1,14) = 0;
          Hf.xelem(2,14) = ((1-r2)*(s+1))/2.0E+0;
          Hf.xelem(0,15) = ((1-r)*(1-s2))/2.0E+0;
          Hf.xelem(1,15) = 0;
          Hf.xelem(2,15) = 0;
          Hf.xelem(0,16) = 0;
          Hf.xelem(1,16) = ((1-r)*(1-s2))/2.0E+0;
          Hf.xelem(2,16) = 0;
          Hf.xelem(0,17) = 0;
          Hf.xelem(1,17) = 0;
          Hf.xelem(2,17) = ((1-r)*(1-s2))/2.0E+0;
          Hf.xelem(0,18) = ((1-r2)*(1-s))/2.0E+0;
          Hf.xelem(1,18) = 0;
          Hf.xelem(2,18) = 0;
          Hf.xelem(0,19) = 0;
          Hf.xelem(1,19) = ((1-r2)*(1-s))/2.0E+0;
          Hf.xelem(2,19) = 0;
          Hf.xelem(0,20) = 0;
          Hf.xelem(1,20) = 0;
          Hf.xelem(2,20) = ((1-r2)*(1-s))/2.0E+0;
          Hf.xelem(0,21) = ((r+1)*(1-s2))/2.0E+0;
          Hf.xelem(1,21) = 0;
          Hf.xelem(2,21) = 0;
          Hf.xelem(0,22) = 0;
          Hf.xelem(1,22) = ((r+1)*(1-s2))/2.0E+0;
          Hf.xelem(2,22) = 0;
          Hf.xelem(0,23) = 0;
          Hf.xelem(1,23) = 0;
          Hf.xelem(2,23) = ((r+1)*(1-s2))/2.0E+0;
     }

     static void VectorInterpMatrixDerR(const ColumnVector& rv, Matrix& dHf_dr) {
          FEM_ASSERT(rv.rows() == 2);

          const double r = rv.xelem(0);
          const double s = rv.xelem(1);
          const double s2 = s * s;

          dHf_dr.xelem(0,0) = (s+1)/4.0E+0-((1-s2)/2.0E+0-r*(s+1))/2.0E+0;
          dHf_dr.xelem(1,0) = 0;
          dHf_dr.xelem(2,0) = 0;
          dHf_dr.xelem(0,1) = 0;
          dHf_dr.xelem(1,1) = (s+1)/4.0E+0-((1-s2)/2.0E+0-r*(s+1))/2.0E+0;
          dHf_dr.xelem(2,1) = 0;
          dHf_dr.xelem(0,2) = 0;
          dHf_dr.xelem(1,2) = 0;
          dHf_dr.xelem(2,2) = (s+1)/4.0E+0-((1-s2)/2.0E+0-r*(s+1))/2.0E+0;
          dHf_dr.xelem(0,3) = (-((-(1-s2)/2.0E+0)-r*(s+1))/2.0E+0)-(s+1)/4.0E+0;
          dHf_dr.xelem(1,3) = 0;
          dHf_dr.xelem(2,3) = 0;
          dHf_dr.xelem(0,4) = 0;
          dHf_dr.xelem(1,4) = (-((-(1-s2)/2.0E+0)-r*(s+1))/2.0E+0)-(s+1)/4.0E+0;
          dHf_dr.xelem(2,4) = 0;
          dHf_dr.xelem(0,5) = 0;
          dHf_dr.xelem(1,5) = 0;
          dHf_dr.xelem(2,5) = (-((-(1-s2)/2.0E+0)-r*(s+1))/2.0E+0)-(s+1)/4.0E+0;
          dHf_dr.xelem(0,6) = (-((-(1-s2)/2.0E+0)-r*(1-s))/2.0E+0)-(1-s)/4.0E+0;
          dHf_dr.xelem(1,6) = 0;
          dHf_dr.xelem(2,6) = 0;
          dHf_dr.xelem(0,7) = 0;
          dHf_dr.xelem(1,7) = (-((-(1-s2)/2.0E+0)-r*(1-s))/2.0E+0)-(1-s)/4.0E+0;
          dHf_dr.xelem(2,7) = 0;
          dHf_dr.xelem(0,8) = 0;
          dHf_dr.xelem(1,8) = 0;
          dHf_dr.xelem(2,8) = (-((-(1-s2)/2.0E+0)-r*(1-s))/2.0E+0)-(1-s)/4.0E+0;
          dHf_dr.xelem(0,9) = (1-s)/4.0E+0-((1-s2)/2.0E+0-r*(1-s))/2.0E+0;
          dHf_dr.xelem(1,9) = 0;
          dHf_dr.xelem(2,9) = 0;
          dHf_dr.xelem(0,10) = 0;
          dHf_dr.xelem(1,10) = (1-s)/4.0E+0-((1-s2)/2.0E+0-r*(1-s))/2.0E+0;
          dHf_dr.xelem(2,10) = 0;
          dHf_dr.xelem(0,11) = 0;
          dHf_dr.xelem(1,11) = 0;
          dHf_dr.xelem(2,11) = (1-s)/4.0E+0-((1-s2)/2.0E+0-r*(1-s))/2.0E+0;
          dHf_dr.xelem(0,12) = -r*(s+1);
          dHf_dr.xelem(1,12) = 0;
          dHf_dr.xelem(2,12) = 0;
          dHf_dr.xelem(0,13) = 0;
          dHf_dr.xelem(1,13) = -r*(s+1);
          dHf_dr.xelem(2,13) = 0;
          dHf_dr.xelem(0,14) = 0;
          dHf_dr.xelem(1,14) = 0;
          dHf_dr.xelem(2,14) = -r*(s+1);
          dHf_dr.xelem(0,15) = -(1-s2)/2.0E+0;
          dHf_dr.xelem(1,15) = 0;
          dHf_dr.xelem(2,15) = 0;
          dHf_dr.xelem(0,16) = 0;
          dHf_dr.xelem(1,16) = -(1-s2)/2.0E+0;
          dHf_dr.xelem(2,16) = 0;
          dHf_dr.xelem(0,17) = 0;
          dHf_dr.xelem(1,17) = 0;
          dHf_dr.xelem(2,17) = -(1-s2)/2.0E+0;
          dHf_dr.xelem(0,18) = -r*(1-s);
          dHf_dr.xelem(1,18) = 0;
          dHf_dr.xelem(2,18) = 0;
          dHf_dr.xelem(0,19) = 0;
          dHf_dr.xelem(1,19) = -r*(1-s);
          dHf_dr.xelem(2,19) = 0;
          dHf_dr.xelem(0,20) = 0;
          dHf_dr.xelem(1,20) = 0;
          dHf_dr.xelem(2,20) = -r*(1-s);
          dHf_dr.xelem(0,21) = (1-s2)/2.0E+0;
          dHf_dr.xelem(1,21) = 0;
          dHf_dr.xelem(2,21) = 0;
          dHf_dr.xelem(0,22) = 0;
          dHf_dr.xelem(1,22) = (1-s2)/2.0E+0;
          dHf_dr.xelem(2,22) = 0;
          dHf_dr.xelem(0,23) = 0;
          dHf_dr.xelem(1,23) = 0;
          dHf_dr.xelem(2,23) = (1-s2)/2.0E+0;
     }

     static void VectorInterpMatrixDerS(const ColumnVector& rv, Matrix& dHf_ds) {
          FEM_ASSERT(rv.rows() == 2);

          const double r = rv.xelem(0);
          const double s = rv.xelem(1);
          const double r2 = r * r;

          dHf_ds.xelem(0,0) = (r+1)/4.0E+0-((1-r2)/2.0E+0-(r+1)*s)/2.0E+0;
          dHf_ds.xelem(1,0) = 0;
          dHf_ds.xelem(2,0) = 0;
          dHf_ds.xelem(0,1) = 0;
          dHf_ds.xelem(1,1) = (r+1)/4.0E+0-((1-r2)/2.0E+0-(r+1)*s)/2.0E+0;
          dHf_ds.xelem(2,1) = 0;
          dHf_ds.xelem(0,2) = 0;
          dHf_ds.xelem(1,2) = 0;
          dHf_ds.xelem(2,2) = (r+1)/4.0E+0-((1-r2)/2.0E+0-(r+1)*s)/2.0E+0;
          dHf_ds.xelem(0,3) = (1-r)/4.0E+0-((1-r2)/2.0E+0-(1-r)*s)/2.0E+0;
          dHf_ds.xelem(1,3) = 0;
          dHf_ds.xelem(2,3) = 0;
          dHf_ds.xelem(0,4) = 0;
          dHf_ds.xelem(1,4) = (1-r)/4.0E+0-((1-r2)/2.0E+0-(1-r)*s)/2.0E+0;
          dHf_ds.xelem(2,4) = 0;
          dHf_ds.xelem(0,5) = 0;
          dHf_ds.xelem(1,5) = 0;
          dHf_ds.xelem(2,5) = (1-r)/4.0E+0-((1-r2)/2.0E+0-(1-r)*s)/2.0E+0;
          dHf_ds.xelem(0,6) = (-((-(1-r)*s)-(1-r2)/2.0E+0)/2.0E+0)-(1-r)/4.0E+0;
          dHf_ds.xelem(1,6) = 0;
          dHf_ds.xelem(2,6) = 0;
          dHf_ds.xelem(0,7) = 0;
          dHf_ds.xelem(1,7) = (-((-(1-r)*s)-(1-r2)/2.0E+0)/2.0E+0)-(1-r)/4.0E+0;
          dHf_ds.xelem(2,7) = 0;
          dHf_ds.xelem(0,8) = 0;
          dHf_ds.xelem(1,8) = 0;
          dHf_ds.xelem(2,8) = (-((-(1-r)*s)-(1-r2)/2.0E+0)/2.0E+0)-(1-r)/4.0E+0;
          dHf_ds.xelem(0,9) = (-((-(r+1)*s)-(1-r2)/2.0E+0)/2.0E+0)-(r+1)/4.0E+0;
          dHf_ds.xelem(1,9) = 0;
          dHf_ds.xelem(2,9) = 0;
          dHf_ds.xelem(0,10) = 0;
          dHf_ds.xelem(1,10) = (-((-(r+1)*s)-(1-r2)/2.0E+0)/2.0E+0)-(r+1)/4.0E+0;
          dHf_ds.xelem(2,10) = 0;
          dHf_ds.xelem(0,11) = 0;
          dHf_ds.xelem(1,11) = 0;
          dHf_ds.xelem(2,11) = (-((-(r+1)*s)-(1-r2)/2.0E+0)/2.0E+0)-(r+1)/4.0E+0;
          dHf_ds.xelem(0,12) = (1-r2)/2.0E+0;
          dHf_ds.xelem(1,12) = 0;
          dHf_ds.xelem(2,12) = 0;
          dHf_ds.xelem(0,13) = 0;
          dHf_ds.xelem(1,13) = (1-r2)/2.0E+0;
          dHf_ds.xelem(2,13) = 0;
          dHf_ds.xelem(0,14) = 0;
          dHf_ds.xelem(1,14) = 0;
          dHf_ds.xelem(2,14) = (1-r2)/2.0E+0;
          dHf_ds.xelem(0,15) = -(1-r)*s;
          dHf_ds.xelem(1,15) = 0;
          dHf_ds.xelem(2,15) = 0;
          dHf_ds.xelem(0,16) = 0;
          dHf_ds.xelem(1,16) = -(1-r)*s;
          dHf_ds.xelem(2,16) = 0;
          dHf_ds.xelem(0,17) = 0;
          dHf_ds.xelem(1,17) = 0;
          dHf_ds.xelem(2,17) = -(1-r)*s;
          dHf_ds.xelem(0,18) = -(1-r2)/2.0E+0;
          dHf_ds.xelem(1,18) = 0;
          dHf_ds.xelem(2,18) = 0;
          dHf_ds.xelem(0,19) = 0;
          dHf_ds.xelem(1,19) = -(1-r2)/2.0E+0;
          dHf_ds.xelem(2,19) = 0;
          dHf_ds.xelem(0,20) = 0;
          dHf_ds.xelem(1,20) = 0;
          dHf_ds.xelem(2,20) = -(1-r2)/2.0E+0;
          dHf_ds.xelem(0,21) = -(r+1)*s;
          dHf_ds.xelem(1,21) = 0;
          dHf_ds.xelem(2,21) = 0;
          dHf_ds.xelem(0,22) = 0;
          dHf_ds.xelem(1,22) = -(r+1)*s;
          dHf_ds.xelem(2,22) = 0;
          dHf_ds.xelem(0,23) = 0;
          dHf_ds.xelem(1,23) = 0;
          dHf_ds.xelem(2,23) = -(r+1)*s;
     }

     static const IntegrationRule& GetIntegrationRule(Element::FemMatrixType eMatType) {
          static constexpr octave_idx_type N = 3;
          static constexpr double r[2][N] = {{0.774596669241483, 0., -0.774596669241483}, {1., 0., -1.}};
          static constexpr double alpha[2][N] = {{0.555555555555556, 0.888888888888889, 0.555555555555556}, {2./3., 2./3., 2./3.}};

          static array<IntegrationRule, 2> rgIntegRule;

          octave_idx_type iIntegRule;

          switch (eMatType) {
          case Element::VEC_LOAD_LUMPED:
               iIntegRule = 1;
               break;
          default:
               iIntegRule = 0;
          }

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

          return rgIntegRule[iIntegRule];
     }
};

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

     static constexpr double EqualityConstr(const ColumnVector&) {
          return 0;
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

     static void ScalarInterpMatrix(const ColumnVector& rv, RowVector& HA) {
          FEM_ASSERT(rv.rows() == 2);

          const double zeta = rv.xelem(0);
          const double eta = rv.xelem(1);

          HA.xelem(0) = (1-2*((-zeta)-eta+1))*(zeta+eta-1);
          HA.xelem(1) = -(1-2*zeta)*zeta;
          HA.xelem(2) = -(1-2*eta)*eta;
          HA.xelem(3) = 4*((-zeta)-eta+1)*zeta;
          HA.xelem(4) = 4*eta*zeta;
          HA.xelem(5) = 4*eta*((-zeta)-eta+1);
     }

     static void VectorInterpMatrix(const ColumnVector& rv, Matrix& Hf) {
          FEM_ASSERT(rv.rows() == 2);

          const double zeta = rv.xelem(0);
          const double eta = rv.xelem(1);

          Hf.xelem(0,0) = (1-2*((-zeta)-eta+1))*(zeta+eta-1);
          Hf.xelem(1,0) = 0;
          Hf.xelem(2,0) = 0;
          Hf.xelem(0,1) = 0;
          Hf.xelem(1,1) = (1-2*((-zeta)-eta+1))*(zeta+eta-1);
          Hf.xelem(2,1) = 0;
          Hf.xelem(0,2) = 0;
          Hf.xelem(1,2) = 0;
          Hf.xelem(2,2) = (1-2*((-zeta)-eta+1))*(zeta+eta-1);
          Hf.xelem(0,3) = -(1-2*zeta)*zeta;
          Hf.xelem(1,3) = 0;
          Hf.xelem(2,3) = 0;
          Hf.xelem(0,4) = 0;
          Hf.xelem(1,4) = -(1-2*zeta)*zeta;
          Hf.xelem(2,4) = 0;
          Hf.xelem(0,5) = 0;
          Hf.xelem(1,5) = 0;
          Hf.xelem(2,5) = -(1-2*zeta)*zeta;
          Hf.xelem(0,6) = -(1-2*eta)*eta;
          Hf.xelem(1,6) = 0;
          Hf.xelem(2,6) = 0;
          Hf.xelem(0,7) = 0;
          Hf.xelem(1,7) = -(1-2*eta)*eta;
          Hf.xelem(2,7) = 0;
          Hf.xelem(0,8) = 0;
          Hf.xelem(1,8) = 0;
          Hf.xelem(2,8) = -(1-2*eta)*eta;
          Hf.xelem(0,9) = 4*((-zeta)-eta+1)*zeta;
          Hf.xelem(1,9) = 0;
          Hf.xelem(2,9) = 0;
          Hf.xelem(0,10) = 0;
          Hf.xelem(1,10) = 4*((-zeta)-eta+1)*zeta;
          Hf.xelem(2,10) = 0;
          Hf.xelem(0,11) = 0;
          Hf.xelem(1,11) = 0;
          Hf.xelem(2,11) = 4*((-zeta)-eta+1)*zeta;
          Hf.xelem(0,12) = 4*eta*zeta;
          Hf.xelem(1,12) = 0;
          Hf.xelem(2,12) = 0;
          Hf.xelem(0,13) = 0;
          Hf.xelem(1,13) = 4*eta*zeta;
          Hf.xelem(2,13) = 0;
          Hf.xelem(0,14) = 0;
          Hf.xelem(1,14) = 0;
          Hf.xelem(2,14) = 4*eta*zeta;
          Hf.xelem(0,15) = 4*eta*((-zeta)-eta+1);
          Hf.xelem(1,15) = 0;
          Hf.xelem(2,15) = 0;
          Hf.xelem(0,16) = 0;
          Hf.xelem(1,16) = 4*eta*((-zeta)-eta+1);
          Hf.xelem(2,16) = 0;
          Hf.xelem(0,17) = 0;
          Hf.xelem(1,17) = 0;
          Hf.xelem(2,17) = 4*eta*((-zeta)-eta+1);
     }

     static void VectorInterpMatrixDerR(const ColumnVector& rv, Matrix& dHf_dr) {
          FEM_ASSERT(rv.rows() == 2);

          const double zeta = rv.xelem(0);
          const double eta = rv.xelem(1);
          
          dHf_dr.xelem(0,0) = 2*(zeta+eta-1)-2*((-zeta)-eta+1)+1;
          dHf_dr.xelem(1,0) = 0;
          dHf_dr.xelem(2,0) = 0;
          dHf_dr.xelem(0,1) = 0;
          dHf_dr.xelem(1,1) = 2*(zeta+eta-1)-2*((-zeta)-eta+1)+1;
          dHf_dr.xelem(2,1) = 0;
          dHf_dr.xelem(0,2) = 0;
          dHf_dr.xelem(1,2) = 0;
          dHf_dr.xelem(2,2) = 2*(zeta+eta-1)-2*((-zeta)-eta+1)+1;
          dHf_dr.xelem(0,3) = 4*zeta-1;
          dHf_dr.xelem(1,3) = 0;
          dHf_dr.xelem(2,3) = 0;
          dHf_dr.xelem(0,4) = 0;
          dHf_dr.xelem(1,4) = 4*zeta-1;
          dHf_dr.xelem(2,4) = 0;
          dHf_dr.xelem(0,5) = 0;
          dHf_dr.xelem(1,5) = 0;
          dHf_dr.xelem(2,5) = 4*zeta-1;
          dHf_dr.xelem(0,6) = 0;
          dHf_dr.xelem(1,6) = 0;
          dHf_dr.xelem(2,6) = 0;
          dHf_dr.xelem(0,7) = 0;
          dHf_dr.xelem(1,7) = 0;
          dHf_dr.xelem(2,7) = 0;
          dHf_dr.xelem(0,8) = 0;
          dHf_dr.xelem(1,8) = 0;
          dHf_dr.xelem(2,8) = 0;
          dHf_dr.xelem(0,9) = 4*((-zeta)-eta+1)-4*zeta;
          dHf_dr.xelem(1,9) = 0;
          dHf_dr.xelem(2,9) = 0;
          dHf_dr.xelem(0,10) = 0;
          dHf_dr.xelem(1,10) = 4*((-zeta)-eta+1)-4*zeta;
          dHf_dr.xelem(2,10) = 0;
          dHf_dr.xelem(0,11) = 0;
          dHf_dr.xelem(1,11) = 0;
          dHf_dr.xelem(2,11) = 4*((-zeta)-eta+1)-4*zeta;
          dHf_dr.xelem(0,12) = 4*eta;
          dHf_dr.xelem(1,12) = 0;
          dHf_dr.xelem(2,12) = 0;
          dHf_dr.xelem(0,13) = 0;
          dHf_dr.xelem(1,13) = 4*eta;
          dHf_dr.xelem(2,13) = 0;
          dHf_dr.xelem(0,14) = 0;
          dHf_dr.xelem(1,14) = 0;
          dHf_dr.xelem(2,14) = 4*eta;
          dHf_dr.xelem(0,15) = -4*eta;
          dHf_dr.xelem(1,15) = 0;
          dHf_dr.xelem(2,15) = 0;
          dHf_dr.xelem(0,16) = 0;
          dHf_dr.xelem(1,16) = -4*eta;
          dHf_dr.xelem(2,16) = 0;
          dHf_dr.xelem(0,17) = 0;
          dHf_dr.xelem(1,17) = 0;
          dHf_dr.xelem(2,17) = -4*eta;
     }

     static void VectorInterpMatrixDerS(const ColumnVector& rv, Matrix& dHf_ds) {
          FEM_ASSERT(rv.rows() == 2);

          const double zeta = rv.xelem(0);
          const double eta = rv.xelem(1);
          
          dHf_ds.xelem(0,0) = 2*(zeta+eta-1)-2*((-zeta)-eta+1)+1;
          dHf_ds.xelem(1,0) = 0;
          dHf_ds.xelem(2,0) = 0;
          dHf_ds.xelem(0,1) = 0;
          dHf_ds.xelem(1,1) = 2*(zeta+eta-1)-2*((-zeta)-eta+1)+1;
          dHf_ds.xelem(2,1) = 0;
          dHf_ds.xelem(0,2) = 0;
          dHf_ds.xelem(1,2) = 0;
          dHf_ds.xelem(2,2) = 2*(zeta+eta-1)-2*((-zeta)-eta+1)+1;
          dHf_ds.xelem(0,3) = 0;
          dHf_ds.xelem(1,3) = 0;
          dHf_ds.xelem(2,3) = 0;
          dHf_ds.xelem(0,4) = 0;
          dHf_ds.xelem(1,4) = 0;
          dHf_ds.xelem(2,4) = 0;
          dHf_ds.xelem(0,5) = 0;
          dHf_ds.xelem(1,5) = 0;
          dHf_ds.xelem(2,5) = 0;
          dHf_ds.xelem(0,6) = 4*eta-1;
          dHf_ds.xelem(1,6) = 0;
          dHf_ds.xelem(2,6) = 0;
          dHf_ds.xelem(0,7) = 0;
          dHf_ds.xelem(1,7) = 4*eta-1;
          dHf_ds.xelem(2,7) = 0;
          dHf_ds.xelem(0,8) = 0;
          dHf_ds.xelem(1,8) = 0;
          dHf_ds.xelem(2,8) = 4*eta-1;
          dHf_ds.xelem(0,9) = -4*zeta;
          dHf_ds.xelem(1,9) = 0;
          dHf_ds.xelem(2,9) = 0;
          dHf_ds.xelem(0,10) = 0;
          dHf_ds.xelem(1,10) = -4*zeta;
          dHf_ds.xelem(2,10) = 0;
          dHf_ds.xelem(0,11) = 0;
          dHf_ds.xelem(1,11) = 0;
          dHf_ds.xelem(2,11) = -4*zeta;
          dHf_ds.xelem(0,12) = 4*zeta;
          dHf_ds.xelem(1,12) = 0;
          dHf_ds.xelem(2,12) = 0;
          dHf_ds.xelem(0,13) = 0;
          dHf_ds.xelem(1,13) = 4*zeta;
          dHf_ds.xelem(2,13) = 0;
          dHf_ds.xelem(0,14) = 0;
          dHf_ds.xelem(1,14) = 0;
          dHf_ds.xelem(2,14) = 4*zeta;
          dHf_ds.xelem(0,15) = 4*((-zeta)-eta+1)-4*eta;
          dHf_ds.xelem(1,15) = 0;
          dHf_ds.xelem(2,15) = 0;
          dHf_ds.xelem(0,16) = 0;
          dHf_ds.xelem(1,16) = 4*((-zeta)-eta+1)-4*eta;
          dHf_ds.xelem(2,16) = 0;
          dHf_ds.xelem(0,17) = 0;
          dHf_ds.xelem(1,17) = 0;
          dHf_ds.xelem(2,17) = 4*((-zeta)-eta+1)-4*eta;
     }

     static const IntegrationRule& GetIntegrationRule(Element::FemMatrixType eMatType) {
          static array<IntegrationRule, 2> rgIntegRule;

          octave_idx_type iIntegRule;
          
          switch (eMatType) {
          case Element::VEC_LOAD_LUMPED:
               iIntegRule = 0;
               break;
          default: 
               iIntegRule = 1;
          }

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
          
          return rgIntegRule[iIntegRule];
     }
};

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

     void SurfaceNormalVector(const Matrix& dHf, ColumnVector& n) const {
          FEM_ASSERT(n.rows() == 3);
          FEM_ASSERT(dHf.rows() == 3);
          FEM_ASSERT(dHf.columns() == nodes.numel() * 3);
          
          for (octave_idx_type i = 0; i < 3; ++i) {
               double ni = 0.;

               for (octave_idx_type j = 0; j < nodes.numel(); ++j) {
                    for (octave_idx_type k = 0; k < 3; ++k) {
                         ni += dHf.xelem(i, j * 3 + k) * X.xelem(k, j);
                    }
               }

               n.xelem(i) = ni;
          }
     }
     
     virtual const IntegrationRule& GetIntegrationRule(FemMatrixType eMatType) const=0;
     virtual void ScalarInterpMatrix(const ColumnVector& r, RowVector& HA) const=0;
     virtual void DisplacementInterpMatrix(const ColumnVector& r, Matrix& Hf) const=0;
     virtual void DisplacementInterpMatrixDerR(const ColumnVector& r, Matrix& dHf_dr) const=0;
     virtual void DisplacementInterpMatrixDerS(const ColumnVector& r, Matrix& dHf_ds) const=0;     
};

class PressureLoad: public SurfaceElement {
public:
     PressureLoad(ElementTypes::TypeId eltype, octave_idx_type id, const Matrix& X, const Material* material, const int32NDArray& nodes, const Matrix& p, octave_idx_type colidx)
          :SurfaceElement(eltype, id, X, material, nodes), p(p), colidx(colidx) {

          FEM_ASSERT(X.rows() == 3);
          FEM_ASSERT(X.columns() == p.columns());
          FEM_ASSERT(X.columns() == nodes.numel());
     }

     PressureLoad(const PressureLoad& oElem)
          :SurfaceElement(oElem.eltype, oElem.id, oElem.X, oElem.material, oElem.nodes), p(oElem.p), colidx(oElem.colidx) {
     }

     virtual void Assemble(MatrixAss& mat, MeshInfo& info, const DofMap& dof, FemMatrixType eMatType) const {
          switch (eMatType) {
          case VEC_LOAD_CONSISTENT:
          case VEC_LOAD_LUMPED:
               break;
          default:
               return;
          }

          const IntegrationRule& oIntegRule = GetIntegrationRule(eMatType);
          const octave_idx_type iNumNodes = nodes.numel();
          const octave_idx_type iNumDof = iGetNumDof();
          const octave_idx_type iNumDir = oIntegRule.iGetNumDirections();
          const octave_idx_type iNumLoads = p.rows();

          ColumnVector rv(iNumDir);
          int32NDArray dofidx(dim_vector(iNumDof, 1), 0);

          for (octave_idx_type inode = 0; inode < iNumNodes; ++inode) {
               for (octave_idx_type idof = 0; idof < 3; ++idof) {
                    dofidx(inode * 3 + idof) = dof.GetNodeDofIndex(nodes(inode).value() - 1, idof);
               }
          }

          RowVector HA(iNumNodes), HA_p(iNumLoads);
          ColumnVector n1(3), n2(3), n_detJA(3), HfT_n_dA(iNumDof);
          Matrix Hf(3, iNumDof), dHf_dr(3, iNumDof), dHf_ds(3, iNumDof), fA(iNumDof, iNumLoads, 0.);

          for (octave_idx_type i = 0; i < oIntegRule.iGetNumEvalPoints(); ++i) {
               const double alpha = oIntegRule.dGetWeight(i);

               for (octave_idx_type j = 0; j < iNumDir; ++j) {
                    rv.xelem(j) = oIntegRule.dGetPosition(i, j);
               }

               DisplacementInterpMatrix(rv, Hf);
               DisplacementInterpMatrixDerR(rv, dHf_dr);
               DisplacementInterpMatrixDerS(rv, dHf_ds);
               ScalarInterpMatrix(rv, HA);

               SurfaceNormalVector(dHf_dr, n1);
               SurfaceNormalVector(dHf_ds, n2);

               n_detJA.xelem(0) = n1.xelem(1) * n2.xelem(2) - n1.xelem(2) * n2.xelem(1);
               n_detJA.xelem(1) = n1.xelem(2) * n2.xelem(0) - n1.xelem(0) * n2.xelem(2);
               n_detJA.xelem(2) = n1.xelem(0) * n2.xelem(1) - n1.xelem(1) * n2.xelem(0);

               double detJA_2 = 0.;

               for (octave_idx_type l = 0; l < 3; ++l) {
                    detJA_2 += n_detJA.xelem(l) * n_detJA.xelem(l);
               }

               info.Add(MeshInfo::JACOBIAN_DET_A, sqrt(detJA_2));

               for (octave_idx_type l = 0; l < iNumDof; ++l) {
                    double HfT_nl_detJA = 0.;

                    for (octave_idx_type m = 0; m < 3; ++m) {
                         HfT_nl_detJA -= Hf.xelem(m, l) * n_detJA.xelem(m);
                    }

                    HfT_n_dA.xelem(l) = HfT_nl_detJA * alpha;
               }

               for (octave_idx_type l = 0; l < iNumLoads; ++l) {
                    double HA_pl = 0.;

                    for (octave_idx_type m = 0; m < iNumNodes; ++m) {
                         HA_pl += HA.xelem(m) * p.xelem(l, m);
                    }

                    HA_p.xelem(l) = HA_pl;
               }

               for (octave_idx_type l = 0; l < iNumLoads; ++l) {
                    for (octave_idx_type m = 0; m < iNumDof; ++m) {
                         fA.xelem(m, l) += HfT_n_dA.xelem(m) * HA_p.xelem(l);
                    }
               }
          }

          for (octave_idx_type j = 0; j < iNumLoads; ++j) {
               for (octave_idx_type i = 0; i < iNumDof; ++i) {
                    mat.Insert(fA.xelem(i, j), dofidx.xelem(i), colidx + j);
               }
          }
     }

     virtual octave_idx_type iGetWorkSpaceSize(FemMatrixType eMatType) const {
          switch (eMatType) {
          case VEC_LOAD_CONSISTENT:
          case VEC_LOAD_LUMPED:
               return iGetNumDof() * p.rows();
          default:
               return 0;
          }
     }

     octave_idx_type iGetNumDof() const {
          return nodes.numel() * 3;
     }

private:
     const Matrix p;
     const octave_idx_type colidx;
};

class ThermalConvection: public SurfaceElement {
public:
     ThermalConvection(ElementTypes::TypeId eltype, octave_idx_type id, const Matrix& X, const Material* material, const int32NDArray& nodes, const Matrix& Thetae, const RowVector& h, octave_idx_type colidx)
          :SurfaceElement(eltype, id, X, material, nodes), Thetae(Thetae), h(h), colidx(colidx) {

          FEM_ASSERT(X.rows() == 3);
          FEM_ASSERT(X.columns() == Thetae.columns());
          FEM_ASSERT(X.columns() == nodes.numel());
     }

     ThermalConvection(const ThermalConvection& oElem)
          :SurfaceElement(oElem.eltype, oElem.id, oElem.X, oElem.material, oElem.nodes), Thetae(oElem.Thetae), h(oElem.h), colidx(oElem.colidx) {
     }

     virtual void Assemble(MatrixAss& mat, MeshInfo& info, const DofMap& dof, const FemMatrixType eMatType) const {
          switch (eMatType) {
          case VEC_LOAD_THERMAL:
               ConvectionLoadVector(mat, info, dof, eMatType);
               break;

          case MAT_THERMAL_COND:
               ConvectionMatrix(mat, info, dof, eMatType);
               break;
               
          default:
               ;
          }          
     }

     void ConvectionMatrix(MatrixAss& mat, MeshInfo& info, const DofMap& dof, const FemMatrixType eMatType) const {
          const IntegrationRule& oIntegRule = GetIntegrationRule(eMatType);
          const octave_idx_type iNumNodes = nodes.numel();
          const octave_idx_type iNumDof = iGetNumDof();
          const octave_idx_type iNumDir = oIntegRule.iGetNumDirections();

          ColumnVector rv(iNumDir);
          int32NDArray dofidx(dim_vector(iNumDof, 1), 0);

          for (octave_idx_type inode = 0; inode < iNumNodes; ++inode) {
               dofidx(inode) = dof.GetNodeDofIndex(nodes(inode).value() - 1, 0);
          }

          RowVector HA(iNumDof);
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

               SurfaceNormalVector(dHf_dr, n1);
               SurfaceNormalVector(dHf_ds, n2);

               double detJA = std::pow(n1.xelem(1) * n2.xelem(2) - n1.xelem(2) * n2.xelem(1), 2)
                    + std::pow(n1.xelem(2) * n2.xelem(0) - n1.xelem(0) * n2.xelem(2), 2)
                    + std::pow(n1.xelem(0) * n2.xelem(1) - n1.xelem(1) * n2.xelem(0), 2);

               if (detJA < 0.) {
                    throw std::runtime_error("Jacobian of surface element is singular");
               }

               detJA = sqrt(detJA);
               
               info.Add(MeshInfo::JACOBIAN_DET_A, detJA);

               double hi = 0.;

               for (octave_idx_type j = 0; j < iNumNodes; ++j) {
                    hi += h.xelem(j) * HA.xelem(j);
               }               
               
               for (octave_idx_type l = 0; l < iNumDof; ++l) {
                    for (octave_idx_type m = l; m < iNumDof; ++m) {
                         KA.xelem(m, l) += hi * HA.xelem(m) * HA.xelem(l) * detJA * alpha;
                    }
               }
          }

          for (octave_idx_type i = 1; i < iNumDof; ++i) {
               for (octave_idx_type j = 0; j < i; ++j) {
                    KA.xelem(j, i) = KA.xelem(i, j);
               }
          }

          for (octave_idx_type j = 0; j < iNumDof; ++j) {
               for (octave_idx_type i = 0; i < iNumDof; ++i) {
                    mat.Insert(KA.xelem(i, j), dofidx.xelem(i), dofidx.xelem(j));
               }
          }
     }
     
     void ConvectionLoadVector(MatrixAss& mat, MeshInfo& info, const DofMap& dof, const FemMatrixType eMatType) const {
          const IntegrationRule& oIntegRule = GetIntegrationRule(eMatType);
          const octave_idx_type iNumNodes = nodes.numel();
          const octave_idx_type iNumDof = iGetNumDof();
          const octave_idx_type iNumDir = oIntegRule.iGetNumDirections();
          const octave_idx_type iNumLoads = Thetae.rows();

          ColumnVector rv(iNumDir);
          int32NDArray dofidx(dim_vector(iNumDof, 1), 0);

          for (octave_idx_type inode = 0; inode < iNumNodes; ++inode) {
               dofidx(inode) = dof.GetNodeDofIndex(nodes(inode).value() - 1, 0);
          }

          RowVector HA(iNumDof), HA_Thetae(iNumLoads);
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

               SurfaceNormalVector(dHf_dr, n1);
               SurfaceNormalVector(dHf_ds, n2);

               double detJA = std::pow(n1.xelem(1) * n2.xelem(2) - n1.xelem(2) * n2.xelem(1), 2)
                    + std::pow(n1.xelem(2) * n2.xelem(0) - n1.xelem(0) * n2.xelem(2), 2)
                    + std::pow(n1.xelem(0) * n2.xelem(1) - n1.xelem(1) * n2.xelem(0), 2);

               if (detJA < 0.) {
                    throw std::runtime_error("Jacobian of surface element is singular");
               }

               detJA = sqrt(detJA);
               
               info.Add(MeshInfo::JACOBIAN_DET_A, detJA);

               double hi = 0.;

               for (octave_idx_type j = 0; j < iNumNodes; ++j) {
                    hi += h.xelem(j) * HA.xelem(j);
               }
               
               for (octave_idx_type l = 0; l < iNumLoads; ++l) {
                    double HA_Thetael = 0.;

                    for (octave_idx_type m = 0; m < iNumNodes; ++m) {
                         HA_Thetael += HA.xelem(m) * Thetae.xelem(l, m);
                    }

                    HA_Thetae.xelem(l) = HA_Thetael;
               }

               for (octave_idx_type l = 0; l < iNumLoads; ++l) {
                    for (octave_idx_type m = 0; m < iNumDof; ++m) {
                         QA.xelem(m, l) += hi * HA.xelem(m) * HA_Thetae.xelem(l) * detJA * alpha;
                    }
               }
          }

          for (octave_idx_type j = 0; j < iNumLoads; ++j) {
               for (octave_idx_type i = 0; i < iNumDof; ++i) {
                    mat.Insert(QA.xelem(i, j), dofidx.xelem(i), colidx + j);
               }
          }
     }

     virtual octave_idx_type iGetWorkSpaceSize(FemMatrixType eMatType) const {
          switch (eMatType) {
          case MAT_THERMAL_COND:
               return iGetNumDof() * iGetNumDof();
               
          case VEC_LOAD_THERMAL:
               return iGetNumDof() * Thetae.rows();
               
          default:
               return 0;
          }
     }

     octave_idx_type iGetNumDof() const {
          return nodes.numel();
     }

private:
     const Matrix Thetae;
     const RowVector h;
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

protected:
     virtual void ScalarInterpMatrix(const ColumnVector& rv, RowVector& HA) const override {
          SHAPE_FUNC::ScalarInterpMatrix(rv, HA);
     }

     virtual void DisplacementInterpMatrix(const ColumnVector& rv, Matrix& Hf) const override {
          SHAPE_FUNC::VectorInterpMatrix(rv, Hf);
     }

     virtual void DisplacementInterpMatrixDerR(const ColumnVector& rv, Matrix& dHf_dr) const override {
          SHAPE_FUNC::VectorInterpMatrixDerR(rv, dHf_dr);
     }

     virtual void DisplacementInterpMatrixDerS(const ColumnVector& rv, Matrix& dHf_ds) const override {
          SHAPE_FUNC::VectorInterpMatrixDerS(rv, dHf_ds);
     }

     virtual const IntegrationRule& GetIntegrationRule(Element::FemMatrixType eMatType) const override {
          return SHAPE_FUNC::GetIntegrationRule(eMatType);
     }
};

typedef SurfaceElementImpl<ShapeIso4, PressureLoad> PressureLoadIso4;
typedef SurfaceElementImpl<ShapeQuad8, PressureLoad> PressureLoadQuad8;
typedef SurfaceElementImpl<ShapeTria6, PressureLoad> PressureLoadTria6;
typedef SurfaceElementImpl<ShapeTria6H, PressureLoad> PressureLoadTria6H;

typedef SurfaceElementImpl<ShapeIso4, ThermalConvection> ThermalConvectionIso4;
typedef SurfaceElementImpl<ShapeQuad8, ThermalConvection> ThermalConvectionQuad8;
typedef SurfaceElementImpl<ShapeTria6, ThermalConvection> ThermalConvectionTria6;
typedef SurfaceElementImpl<ShapeTria6H, ThermalConvection> ThermalConvectionTria6H;

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

     virtual void Assemble(MatrixAss& mat, MeshInfo& info, const DofMap& dof, FemMatrixType eMatType) const {
          switch (eMatType) {
          case VEC_LOAD_CONSISTENT:
          case VEC_LOAD_LUMPED:
               break;
          default:
               return;
          }

          for (octave_idx_type j = 0; j < loads.columns(); ++j) {
               for (octave_idx_type i = 0; i < loads.rows(); ++i) {
                    const octave_idx_type inode = nodes.xelem(i).value() - 1;

                    mat.Insert(loads.xelem(i, j), dof.GetNodeDofIndex(inode, j), colidx);
               }
          }
     }

     virtual octave_idx_type iGetWorkSpaceSize(FemMatrixType eMatType) const {
          switch (eMatType) {
          case VEC_LOAD_CONSISTENT:
          case VEC_LOAD_LUMPED:
               return loads.rows() * loads.columns();
          default:
               return 0;
          }
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
     virtual void Assemble(MatrixAss& oMatAss, MeshInfo& info, const DofMap& oDof, Element::FemMatrixType eMatType) const=0;
     virtual void PostProcElem(NDArray& mat, Element::FemMatrixType eMatType, const NDArray& U) const=0;
     virtual double dGetMass() const=0;
     virtual bool bNeedMatrixInfo(Element::FemMatrixType eMatType) const=0;

     template <typename ElementType, typename... Args>
     inline void Insert(octave_idx_type id, const Matrix& X, const Material* material, const int32NDArray& nodes, const Args&... args);

     ElementTypes::TypeId GetElementType() const { return eltype; }
     virtual void Extract(octave_idx_type& idx, octave_map& sElem) const=0;
     virtual octave_idx_type iGetNumElem() const=0;

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
          rgElements.reserve(elements.rows());

          Matrix X_e(inumcoord, elements.columns());
          int32NDArray nodes_e(dim_vector(elements.columns(), 1));

          for (octave_idx_type i = 0; i < elements.rows(); ++i) {
               X_e.make_unique();
               nodes_e.make_unique();

               for (octave_idx_type j = 0; j < elements.columns(); ++j) {
                    nodes_e.xelem(j) = elements.xelem(i, j);
                    for (octave_idx_type k = 0; k < X_e.rows(); ++k) {
                         X_e.xelem(k, j) = nodes.xelem(nodes_e.xelem(j).value() - 1, k);
                    }
               }

               const Material* material = nullptr; // Some elements like RBE3 do not need a material

               if (materials(i).value() > 0) {
                    material = &rgMaterials[materials.xelem(i).value() - 1];
               }

               rgElements.emplace_back(eltype, i + 1, X_e, material, nodes_e, args...);
          }
     }

     template <typename... Args>
     void Insert(octave_idx_type id, const Matrix& X, const Material* material, const int32NDArray& nodes, const Args&... args) {
          rgElements.emplace_back(eltype, id, X, material, nodes, args...);
     }

     octave_idx_type iGetWorkSpaceSize(Element::FemMatrixType eMatType) const {
          octave_idx_type iWorkSpace = 0;

          for (auto i = rgElements.begin(); i != rgElements.end(); ++i) {
               iWorkSpace += i->ElementType::iGetWorkSpaceSize(eMatType);
          }

          return iWorkSpace;
     }

     void Assemble(MatrixAss& oMatAss, MeshInfo& oMeshInfo, const DofMap& oDof, Element::FemMatrixType eMatType) const {
          for (auto i = rgElements.begin(); i != rgElements.end(); ++i) {
               i->ElementType::Assemble(oMatAss, oMeshInfo, oDof, eMatType);

               OCTAVE_QUIT;
          }
     }

     void PostProcElem(NDArray& mat, Element::FemMatrixType eMatType, const NDArray& U) const {
          for (auto i = rgElements.begin(); i != rgElements.end(); ++i) {
               i->ElementType::PostProcElem(mat, eMatType, U);

               OCTAVE_QUIT;
          }
     }

     double dGetMass() const {
          double dm = 0.;

          for (auto i = rgElements.begin(); i != rgElements.end(); ++i) {
               dm += i->ElementType::dGetMass();

               OCTAVE_QUIT;
          }

          return dm;
     }

     virtual bool bNeedMatrixInfo(Element::FemMatrixType eMatType) const {
          return rgElements.empty() ? false : rgElements.front().bNeedMatrixInfo(eMatType);
     }

     void Reserve(octave_idx_type iNumElem) {
          rgElements.reserve(iNumElem);
     }

     virtual octave_idx_type iGetNumElem() const {
          return rgElements.size();
     }

     virtual void Extract(octave_idx_type& idx, octave_map& sElem) const {
          for (const auto& oElem: rgElements) {
               FEM_ASSERT(sElem.numel() > idx);
               oElem.ElementType::Extract(idx, sElem);
          }
     }
private:
     vector<ElementType> rgElements;
};

template <typename ElementType, typename... Args>
void ElementBlockBase::Insert(octave_idx_type id, const Matrix& X, const Material* material, const int32NDArray& nodes, const Args&... args) {
     typedef ElementBlock<ElementType> ElemBlockType;

     FEM_ASSERT(dynamic_cast<ElemBlockType*>(this) == static_cast<ElemBlockType*>(this));

     static_cast<ElemBlockType*>(this)->Insert(id, X, material, nodes, args...);
}

template <typename PressureElemType>
void InsertPressureElem(ElementTypes::TypeId eltype, const Matrix& nodes, const octave_map& load_case, const char* pszElemName, octave_idx_type iNumNodesElem, vector<std::unique_ptr<ElementBlockBase> >& rgElemBlocks) {
     const auto iter_pressure = load_case.seek("pressure");

     if (iter_pressure != load_case.end()) {
          const Cell cell_pressure = load_case.contents(iter_pressure);

          FEM_ASSERT(cell_pressure.numel() == load_case.numel());

          std::unique_ptr<ElementBlock<PressureElemType> > pElem(nullptr);

          for (octave_idx_type j = 0; j < 2; ++j) {
               octave_idx_type iNumElements = 0;

               for (octave_idx_type i = 0; i < cell_pressure.numel(); ++i) {
                    if (cell_pressure(i).isstruct()) {
                         if (!(cell_pressure(i).numel() == 1)) {
                              throw std::runtime_error("pressure must be a scalar struct");
                         }

                         const octave_scalar_map pressure = cell_pressure.xelem(i).scalar_map_value();

                         const auto iter_elem_type = pressure.seek(pszElemName);

                         if (iter_elem_type == pressure.end()) {
                              continue;
                         }

                         const octave_value ov_elem_type = pressure.contents(iter_elem_type);

                         if (!(ov_elem_type.isstruct() && ov_elem_type.numel() == 1)) {
                              throw std::runtime_error("invalid entry in load_case.pressure");
                         }

                         const octave_scalar_map elem_type = ov_elem_type.scalar_map_value();

                         const auto iter_elements = elem_type.seek("elements");

                         if (iter_elements == elem_type.end()) {
                              throw std::runtime_error("field \"elements\" not found in struct pressure");
                         }

                         const auto iter_p = elem_type.seek("p");

                         if (iter_p == elem_type.end()) {
                              throw std::runtime_error("field \"p\" not found in struct pressure");
                         }

                         const octave_value ov_elements = elem_type.contents(iter_elements);

                         if (!(ov_elements.is_matrix_type() && ov_elements.OV_ISINTEGER())) {
                              throw std::runtime_error("pressure.elements must be an integer array");
                         }

                         const int32NDArray elements = ov_elements.int32_array_value();

                         if (elements.columns() != iNumNodesElem) {
                              throw std::runtime_error("pressure.elements number of columns do not match");
                         }

                         const octave_value ov_p = elem_type.contents(iter_p);

                         if (!(ov_p.is_matrix_type() && ov_p.OV_ISREAL())) {
                              throw std::runtime_error("pressure.p must be a real matrix");
                         }

                         const Matrix p = ov_p.matrix_value();

                         if (p.columns() != elements.columns() || p.rows() != elements.rows()) {
                              throw std::runtime_error("pressure.p must have the same shape like pressure.elements");
                         }

                         switch (j) {
                         case 0:
                              iNumElements += elements.rows();
                              break;

                         case 1: {
                              Matrix X(3, iNumNodesElem);

                              for (octave_idx_type k = 0; k < elements.rows(); ++k) {
                                   X.make_unique();

                                   for (octave_idx_type l = 0; l < X.columns(); ++l) {
                                        for (octave_idx_type m = 0; m < X.rows(); ++m) {
                                             octave_idx_type inode = elements.xelem(k, l).value() - 1;

                                             if (inode < 0 || inode >= nodes.rows()) {
                                                  throw std::runtime_error("node index out of range in pressure.elements");
                                             }

                                             X.xelem(m, l) = nodes.xelem(inode, m);
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
void InsertThermalConvElem(ElementTypes::TypeId eltype, const Matrix& nodes, const octave_scalar_map& elements, const octave_map& load_case, const char* pszElemName, octave_idx_type iNumNodesElem, vector<std::unique_ptr<ElementBlockBase> >& rgElemBlocks) {
     const auto iter_convection = elements.seek("convection");

     if (iter_convection == elements.end()) {
          return;
     }
     
     const octave_value ov_convection = elements.contents(iter_convection);

     if (!(ov_convection.isstruct() && ov_convection.numel() == 1)) {
          throw std::runtime_error("mesh.elements.convection must be a scalar struct");
     }

     const octave_scalar_map m_convection = ov_convection.scalar_map_value();
     
     const auto iter_elem_type = m_convection.seek(pszElemName);

     if (iter_elem_type == m_convection.end()) {
          return;
     }

     const octave_value ov_elem_type = m_convection.contents(iter_elem_type);

     if (!(ov_elem_type.isstruct() && ov_elem_type.numel() == 1)) {
          throw std::runtime_error("mesh.elements.convection."s + pszElemName + " must be a scalar struct");
     }

     const octave_scalar_map m_elem_type = ov_elem_type.scalar_map_value();

     const auto iter_elnodes = m_elem_type.seek("nodes");

     if (iter_elnodes == m_elem_type.end()) {
          throw std::runtime_error("missing field mesh.elements.convection."s + pszElemName + ".nodes");
     }

     const octave_value ov_elnodes = m_elem_type.contents(iter_elnodes);

     if (!(ov_elnodes.is_matrix_type() && ov_elnodes.isinteger() && ov_elnodes.columns() == iNumNodesElem)) {
          throw std::runtime_error("mesh.elements.convection."s + pszElemName + ".nodes must be an integer matrix");          
     }

     const int32NDArray elnodes = ov_elnodes.int32_array_value();

     const auto iter_h = m_elem_type.seek("h");

     if (iter_h == m_elem_type.end()) {
          throw std::runtime_error("missing field mesh.elements.convection."s + pszElemName + ".h");
     }

     const octave_value ov_h = m_elem_type.contents(iter_h);

     if (!(ov_h.is_matrix_type() && ov_h.isreal() && ov_h.rows() == elnodes.rows() && ov_h.columns() == elnodes.columns())) {
          throw std::runtime_error("mesh.elements.convection."s + pszElemName
                                   + ".h must be a real matrix with the same size like mesh.elements.convection."
                                   + pszElemName + ".nodes");
     }

     const Matrix h = ov_h.matrix_value();

     NDArray theta(dim_vector(load_case.numel(), elnodes.columns(), elnodes.rows()), 0.);

     const auto iter_conv_load = load_case.seek("convection");

     if (iter_conv_load != load_case.end()) {
          const Cell cell_conv_load = load_case.contents(iter_conv_load);

          FEM_ASSERT(cell_conv_load.numel() == theta.rows());
          
          for (octave_idx_type k = 0; k < theta.rows(); ++k) {
               const octave_value ov_conv_load = cell_conv_load.xelem(k);
               
               if (!(ov_conv_load.isstruct() && ov_conv_load.numel() == 1)) {
                    throw std::runtime_error("load_case.convection must be a scalar struct");
               }

               const octave_scalar_map m_conv_load = ov_conv_load.scalar_map_value();

               const auto iter_conv_load_elem = m_conv_load.seek(pszElemName);

               if (iter_conv_load_elem == m_conv_load.end()) {
                    continue;
               }

               const octave_value ov_conv_load_elem = m_conv_load.contents(iter_conv_load_elem);

               if (!(ov_conv_load_elem.isstruct() && ov_conv_load_elem.numel() == 1)) {
                    throw std::runtime_error("load_case.convection."s + pszElemName + " must be a scalar struct");
               }

               const octave_scalar_map m_conv_load_elem = ov_conv_load_elem.scalar_map_value();

               const auto iter_theta = m_conv_load_elem.seek("theta");

               if (iter_theta == m_conv_load_elem.end()) {
                    throw std::runtime_error("missing field load_case.convection."s + pszElemName + ".theta");
               }

               const octave_value ov_thetak = m_conv_load_elem.contents(iter_theta);

               if (!(ov_thetak.is_matrix_type() && ov_thetak.isreal() && ov_thetak.rows() == elnodes.rows() && ov_thetak.columns() == elnodes.columns())) {
                    throw std::runtime_error("load_case.convection."s + pszElemName + ".theta must be a real matrix with the same dimensions like mesh.elements.convection."s + pszElemName + ".nodes");
               }

               const Matrix thetak = ov_thetak.matrix_value();

               for (octave_idx_type j = 0; j < elnodes.rows(); ++j) {
                    for (octave_idx_type i = 0; i < elnodes.columns(); ++i) {
                         theta.xelem(k, i, j) = thetak.xelem(j, i);
                    }
               }
          }
     }
     
     NDArray X(dim_vector(3, iNumNodesElem, elnodes.rows()));     

     for (octave_idx_type k = 0; k < elnodes.rows(); ++k) {
          for (octave_idx_type l = 0; l < elnodes.columns(); ++l) {
               octave_idx_type inode = elnodes.xelem(k, l).value() - 1;

               if (inode < 0 || inode >= nodes.rows()) {
                    throw std::runtime_error("node index out of range in mesh.elements.convection."s + pszElemName + ".nodes");
               }
               
               for (octave_idx_type m = 0; m < X.rows(); ++m) {
                    X.xelem(m, l, k) = nodes.xelem(inode, m);
               }
          }
     }
     
     std::unique_ptr<ElementBlock<ConvectionElemType> > pElem{new ElementBlock<ConvectionElemType>(eltype)};

     pElem->Reserve(elnodes.rows());

     for (octave_idx_type k = 0; k < elnodes.rows(); ++k) {
          const Matrix Xk = X.linear_slice(X.rows() * X.columns() * k, X.rows() * X.columns() * (k + 1)).reshape(dim_vector(X.rows(), X.columns()));
          const Matrix thetak = theta.linear_slice(theta.rows() * theta.columns() * k, theta.rows() * theta.columns() * (k + 1)).reshape(dim_vector(theta.rows(), theta.columns()));
          
          pElem->Insert(k,
                        Xk,
                        nullptr,
                        elnodes.index(idx_vector::make_range(k, 1, 1),
                                      idx_vector::make_range(0, 1, elnodes.columns())),
                        thetak,
                        h.row(k),
                        1);
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
     struct SelectElemPerSlaveNode<ShapeTria6> {
          static constexpr octave_idx_type iNumElemPerSlaveNode = 13;
     };

     template <>
     struct SelectElemPerSlaveNode<ShapeTria6H> {
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

     typedef std::unique_ptr<ElementBlock<ElemJoint>> ElemBlockPtr;

     static ConstraintType GetConstraintType(int constr) {
          switch (constr) {
          case CT_FIXED:
          case CT_SLIDING:
               return static_cast<ConstraintType>(constr);
          default:
               throw std::runtime_error("invalid value for elements.sfncon{4|6|8}.constraint");
          }
     }

     static ConstraintType GetConstraintType(const Cell& ov, octave_idx_type j) {
          if (!ov.numel()) {
               return CT_FIXED;
          }

          const int constr = ov.xelem(j).int_value();

          if (error_state) {
               throw std::runtime_error("elements.sfncon{4|6|8}.constraint must be a scalar value");
          }

          return GetConstraintType(constr);
     }

     static octave_idx_type iGetNumDof(ConstraintType eType) {
          return eType == CT_FIXED ? 3 : 1;
     }

     static octave_idx_type iGetNumDof(const Cell& ov, octave_idx_type j) {
          return iGetNumDof(GetConstraintType(ov, j));
     }

     static void BuildJoints(const Matrix& nodes,
                             const octave_scalar_map& elements,
                             const array<int32NDArray, DofMap::ELEM_TYPE_COUNT>& edof,
                             array<octave_idx_type, DofMap::ELEM_TYPE_COUNT>& dofelemid,
                             const ElementTypes::TypeInfo& oElemType,
                             const unsigned uConstraintFlags,
                             vector<std::unique_ptr<ElementBlockBase> >& rgElemBlocks);
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
     static ElemBlockPtr
     BuildJoints(octave_idx_type& id,
                 const Matrix& X,
                 const int32NDArray& nidxmaster,
                 const int32NDArray& nidxslave,
                 const ColumnVector& maxdist,
                 const ConstraintType eType,
                 const unsigned uConstraintFlags) {
          ElemBlockPtr pElemBlock(new ElementBlock<ElemJoint>(ElementTypes::ELEM_JOINT, nidxslave.numel()));

          FEM_ASSERT(X.columns() >= iNumDofNode);
          FEM_ASSERT(X.columns() == 6);
          FEM_ASSERT(nidxmaster.columns() == iNumNodesElem);
          FEM_ASSERT(nidxslave.rows() == maxdist.rows());

          ColumnVector Xs(iNumDofNode);

          vector<ElemIndexVector> eidxmaster(nidxslave.numel());

          for (octave_idx_type k = 0; k < nidxslave.numel(); ++k) {
               for (octave_idx_type l = 0; l < Xs.rows(); ++l) {
                    Xs.xelem(l) = X.xelem(nidxslave.xelem(k).value() - 1, l);
               }

               for (octave_idx_type i = 0; i < nidxmaster.rows(); ++i) {
                    double dXmin = std::numeric_limits<double>::max();

                    for (octave_idx_type j = 0; j < nidxmaster.columns(); ++j) {
                         double dX = 0;

                         for (octave_idx_type l = 0; l < Xs.rows(); ++l) {
                              dX += std::pow(X.xelem(nidxmaster.xelem(i, j).value() - 1, l) - Xs.xelem(l), 2);
                         }

                         dXmin = std::min(dX, dXmin);
                    }

                    eidxmaster[k].insert(ElemIndexRecord(dXmin, i));
               }

               OCTAVE_QUIT;
          }

          ColumnVector Xm(iNumDofNode * iNumNodesElem);
          Matrix Xe(6, nidxmaster.columns() + 1);
          ColumnVector rv(iNumDir), rvopt(iNumDir, 0.);
          ColumnVector Xi(iNumDofNode), Xiopt(iNumDofNode, 0);
          Matrix Hf(iNumDofNode, iNumDofNode * iNumNodesElem);
          Matrix dHf_dr(iNumDofNode, iNumDofNode * iNumNodesElem);
          Matrix dHf_ds(iNumDofNode, iNumDofNode * iNumNodesElem);
          ColumnVector n1(3), n2(3);
          RowVector n(3);

          Matrix C(iNumDofNode, 6 * (nidxmaster.columns() + 1));
          Matrix U(C.rows(), 0);
          RowVector nC(C.columns());
          Matrix nU(1, 0);
          int32NDArray enodes(dim_vector(nidxmaster.columns() + 1, 1));

          for (size_t i = 0; i < eidxmaster.size(); ++i) {
               for (octave_idx_type j = 0; j < Xs.rows(); ++j) {
                    Xs.xelem(j) = X.xelem(nidxslave.xelem(i).value() - 1, j);
               }

               double fopt = std::numeric_limits<double>::max();
               nlopt_result rcopt = NLOPT_FAILURE;
               octave_idx_type lopt = -1;
               rvopt.fill(0.);

               for (octave_idx_type l = 0; l < eidxmaster[i].size(); ++l) {
                    for (octave_idx_type j = 0; j < nidxmaster.columns(); ++j) {
                         for (octave_idx_type k = 0; k < iNumDofNode; ++k) {
                              Xm.xelem(j * iNumDofNode + k) = X.xelem(nidxmaster.xelem(eidxmaster[i][l].eidx, j).value() - 1, k);
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

               if (rcopt < 0 || fopt > maxdist(i) || lopt < 0) {
                    if (uConstraintFlags & CF_IGNORE_NODES_OUT_OF_RANGE) {
                         continue;
                    }

                    std::ostringstream os;

                    os << "nlopt failed to project slave node #" << i + 1 << " (" << nidxslave(i).value()
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

                    for (octave_idx_type k = 0; k < iNumDofNode; ++k) {
                         os << X.xelem(nidxslave.xelem(i).value() - 1, k) << ' ';
                    }

                    os << std::endl;

                    os << "Xi=" << Xiopt.transpose() << std::endl;
                    os << "Xm=" << std::endl;

                    if (eidx >= 0) {
                         for (octave_idx_type k = 0; k < iNumNodesElem; ++k) {
                              for (octave_idx_type j = 0; j < iNumDofNode; ++j) {
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

               SHAPE_FUNC::VectorInterpMatrix(rvopt, Hf);

               C.make_unique();
               C.fill(0.);

               for (octave_idx_type j = 0; j < 6; ++j) {
                    for (octave_idx_type k = 0; k < C.rows(); ++k) {
                         C.xelem(k, j) = (k == j);
                    }
               }

               for (octave_idx_type j = 0; j < nidxmaster.columns(); ++j) {
                    for (octave_idx_type k = 0; k < iNumDofNode; ++k) {
                         for (octave_idx_type l = 0; l < C.rows(); ++l) {
                              C.xelem(l, (j + 1) * 6 + k) = -Hf.xelem(l, j * iNumDofNode + k);
                         }
                    }
               }

               enodes.make_unique();
               enodes.xelem(0) = nidxslave.xelem(i);

               for (octave_idx_type j = 0; j < nidxmaster.columns(); ++j) {
                    enodes.xelem(j + 1) = nidxmaster.xelem(eidxmaster[i][lopt].eidx, j).value();
               }

               for (octave_idx_type j = 0; j < X.columns(); ++j) {
                    Xe.xelem(j, 0) = X.xelem(nidxslave.xelem(i).value() - 1, j);
               }

               for (octave_idx_type j = 0; j < nidxmaster.columns(); ++j) {
                    for (octave_idx_type k = 0; k < X.columns(); ++k) {
                         Xe.xelem(k, j + 1) = X.xelem(nidxmaster.xelem(eidxmaster[i][lopt].eidx, j).value() - 1, k);
                    }
               }

               switch (eType) {
               case CT_SLIDING: {
                    SHAPE_FUNC::VectorInterpMatrixDerR(rvopt, dHf_dr);
                    SHAPE_FUNC::VectorInterpMatrixDerS(rvopt, dHf_ds);
                    SurfaceNormalVector(Xe, dHf_dr, n1);
                    SurfaceNormalVector(Xe, dHf_ds, n2);

                    n.xelem(0) = n1.xelem(1) * n2.xelem(2) - n1.xelem(2) * n2.xelem(1);
                    n.xelem(1) = n1.xelem(2) * n2.xelem(0) - n1.xelem(0) * n2.xelem(2);
                    n.xelem(2) = n1.xelem(0) * n2.xelem(1) - n1.xelem(1) * n2.xelem(0);

                    double norm_n = 0.;

                    for (octave_idx_type l = 0; l < 3; ++l) {
                         norm_n += n.xelem(l) * n(l);
                    }

                    if (norm_n <= 0) {
                         std::ostringstream os;

                         os << "Jacobian of sfncon"
                            << SHAPE_FUNC::iGetNumNodes()
                            << " element " << id << " is singular"
                            << std::ends;

                         throw std::runtime_error(os.str());
                    }

                    norm_n = sqrt(norm_n);

                    for (octave_idx_type l = 0; l < 3; ++l) {
                         n.xelem(l) /= norm_n;
                    }

                    nC.make_unique();
                    nC.fill(0.);

                    for (octave_idx_type j = 0; j < C.columns(); ++j) {
                         for (octave_idx_type k = 0; k < C.rows(); ++k) {
                              nC.xelem(j) += n.xelem(k) * C.xelem(k, j);
                         }
                    }

                    FEM_TRACE("n=" << n << std::endl);
                    FEM_TRACE("C=" << C << std::endl);
                    FEM_TRACE("n.'*C=" << nC << std::endl);

                    pElemBlock->Insert(++id, Xe, nullptr, enodes, nC, nU);
               } break;
               default:
                    pElemBlock->Insert(++id, Xe, nullptr, enodes, C, U);
                    break;
               }
          }

          return pElemBlock;
     }

private:
     static void SurfaceNormalVector(const Matrix& X, const Matrix& dHf, ColumnVector& n) {
          for (octave_idx_type i = 0; i < 3; ++i) {
               double ni = 0.;

               for (octave_idx_type j = 0; j < X.columns() - 1; ++j) {
                    for (octave_idx_type k = 0; k < 3; ++k) {
                         ni += dHf.xelem(i, j * 3 + k) * X.xelem(k, j + 1);
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
               nlopt_add_equality_constraint(oFuncData.opt, &SurfToNodeConstr::EqualityConstr, &oFuncData, dTolX);
          }

          Matrix dX(iNumDofNode, 2);

          for (octave_idx_type i = 0; i < iNumDofNode; ++i) {
               dX.xelem(i, 0) = std::numeric_limits<double>::max();
               dX.xelem(i, 1) = -dX.xelem(i, 0);

               for (octave_idx_type j = 0; j < iNumNodesElem; ++j) {
                    dX.xelem(i, 0) = std::min(dX.xelem(i, 0), Xm.xelem(j * iNumDofNode + i));
                    dX.xelem(i, 1) = std::max(dX.xelem(i, 1), Xm.xelem(j * iNumDofNode + i));
               }
          }

          double dTolF = 0;

          for (octave_idx_type i = 0; i < iNumDofNode; ++i) {
               dTolF = std::max(dTolF, dTolX * (dX.xelem(i, 1) - dX.xelem(i, 0)));
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

     static constexpr octave_idx_type iNumDofNode = SHAPE_FUNC::iGetNumDofNode();
     static constexpr octave_idx_type iNumDir = SHAPE_FUNC::iGetNumDirections();
     static constexpr octave_idx_type iNumNodesElem = SHAPE_FUNC::iGetNumNodes();

     SurfToNodeConstr(const ColumnVector& Xm, const ColumnVector& Xs)
          :Xm(Xm),
           Xs(Xs) {

          FEM_ASSERT(Xm.rows() == iNumNodesElem * Xs.rows());
          FEM_ASSERT(Xs.rows() == iNumDofNode);
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

     double Objective(const ColumnVector& rv, ColumnVector& f, Matrix& Hf) const {
          FEM_ASSERT(f.rows() == iNumDofNode);
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

          FEM_TRACE("f=" << sqrt(ftot) << std::endl);

          return ftot;
     }

     static double EqualityConstr(const ColumnVector& rv) {
          FEM_ASSERT(SHAPE_FUNC::iGetNumEqualityConstr());

          double f = SHAPE_FUNC::EqualityConstr(rv);

          OCTAVE_QUIT;

          return f;
     }

     void Position(const ColumnVector& rv, ColumnVector& Xi, Matrix& Hf) const {
          FEM_ASSERT(rv.rows() == iNumDir);
          FEM_ASSERT(Xi.rows() == iNumDofNode);
          FEM_ASSERT(Hf.rows() == iNumDofNode);
          FEM_ASSERT(Hf.columns() == Xi.rows() * iNumNodesElem);
          FEM_ASSERT(Xm.rows() == Hf.columns());

          SHAPE_FUNC::VectorInterpMatrix(rv, Hf);

          for (octave_idx_type i = 0; i < Xi.rows(); ++i) {
               Xi.xelem(i) = 0;
          }

          for (octave_idx_type j = 0; j < Hf.columns(); ++j) {
               for (octave_idx_type i = 0; i < Hf.rows(); ++i) {
                    Xi.xelem(i) += Hf.xelem(i, j) * Xm.xelem(j);
               }
          }
     }

private:
     const ColumnVector Xm, Xs;

     struct FuncData {
          FuncData(const ColumnVector& Xm, const ColumnVector& Xs)
               :oSNCO(Xm, Xs),
                opt(nlopt_create(NLOPT_LD_SLSQP, iNumDir)),
                Hf(iNumDofNode, iNumNodesElem * iNumDofNode),
                rv(iNumDir, 0.),
                f(iNumDofNode) {
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
                                       vector<std::unique_ptr<ElementBlockBase> >& rgElemBlocks) {
     const auto iter_elem = elements.seek(oElemType.name);

     if (iter_elem == elements.end()) {
          return;
     }

     const octave_map s_elem(elements.contents(iter_elem).map_value());

     if (error_state) {
          throw std::runtime_error("elements.sfncon{4|6|8} must be a struct array");
     }

     const auto iter_nidxmaster = s_elem.seek("master");

     if (iter_nidxmaster == s_elem.end()) {
          throw std::runtime_error("elements.sfncon{4|6|8}.master not defined");
     }

     const Cell ov_nidxmaster = s_elem.contents(iter_nidxmaster);

     const auto iter_nidxslave = s_elem.seek("slave");

     if (iter_nidxslave == s_elem.end()) {
          throw std::runtime_error("elements.sfncon{4|6|8}.slave not defined");
     }

     const Cell ov_nidxslave = s_elem.contents(iter_nidxslave);

     const auto iter_maxdist = s_elem.seek("maxdist");

     if (iter_maxdist == s_elem.end()) {
          throw std::runtime_error("elements.sfncon{4|6|8}.maxdist not defined");
     }

     const Cell ov_maxdist = s_elem.contents(iter_maxdist);

     const auto iter_constr = s_elem.seek("constraint");

     Cell ov_constr;

     if (iter_constr != s_elem.end()) {
          ov_constr = s_elem.contents(iter_constr);
     }

     FEM_ASSERT(ov_nidxslave.numel() == s_elem.numel());
     FEM_ASSERT(ov_nidxmaster.numel() == s_elem.numel());
     FEM_ASSERT(ov_maxdist.numel() == s_elem.numel());
     FEM_ASSERT(ov_constr.numel() == 0 || ov_constr.numel() == s_elem.numel());

     for (octave_idx_type l = 0; l < s_elem.numel(); ++l) {
          const int32NDArray nidxmaster = ov_nidxmaster(l).int32_array_value();

          if (error_state) {
               throw std::runtime_error("elements.sfncon{4|6|8}.master must be an integer array");
          }

          const int32NDArray nidxslave = ov_nidxslave(l).int32_array_value();

          if (error_state) {
               throw std::runtime_error("elements.sfncon{4|6|8}.slave must be an integer array");
          }

          ColumnVector maxdist = ov_maxdist(l).column_vector_value();

          if (error_state) {
               throw std::runtime_error("elements.sfncon{4|6|8}.maxdist must be a column vector");
          }

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
          default:
               FEM_ASSERT(false);
          }

          if (nidxmaster.ndims() != 2 || nidxmaster.rows() < 1 || nidxmaster.columns() != iNumNodesElem) {
               throw std::runtime_error("elements.sfncon{4|6|8}.master must be an nx{4|6|8} array");
          }

          if (nidxslave.ndims() != 2 || nidxslave.rows() < 1 || nidxslave.columns() != 1) {
               throw std::runtime_error("elements.sfncon{4|6|8}.slave must be an nx1 array");
          }

          if (maxdist.rows() == 1 && nidxslave.rows() > 1) {
               const double maxdistval = maxdist(0);
               maxdist.resize(nidxslave.rows(), maxdistval);
          }

          if (maxdist.rows() != nidxslave.rows()) {
               throw std::runtime_error("elements.sfncon{4|6|8}.maxdist must have "
                                        "the same dimensions like elements.sfncon{4|6|8}.slave");
          }

          for (octave_idx_type i = 0; i < nidxmaster.rows(); ++i) {
               for (octave_idx_type j = 0; j < nidxmaster.columns(); ++j) {
                    if (nidxmaster(i, j).value() < 1 || nidxmaster(i, j).value() > nodes.rows()) {
                         throw std::runtime_error("elements.sfncon{4|6|8}.master: node index out of range");
                    }
               }
          }

          for (octave_idx_type i = 0; i < nidxslave.rows(); ++i) {
               if (nidxslave(i).value() < 1 || nidxslave(i).value() > nodes.rows()) {
                    throw std::runtime_error("elements.sfncon{4|6|8}.slave: node index out of range");
               }
          }

          std::unique_ptr<ElementBlock<ElemJoint> > pElem;

          switch (oElemType.type) {
          case ElementTypes::ELEM_SFNCON4:
               pElem = SurfToNodeConstr<ShapeIso4>::BuildJoints(dofelemid[oElemType.dof_type],
                                                                nodes,
                                                                nidxmaster,
                                                                nidxslave,
                                                                maxdist,
                                                                eConstrType,
                                                                uConstraintFlags);
               break;
          case ElementTypes::ELEM_SFNCON6:
               pElem = SurfToNodeConstr<ShapeTria6>::BuildJoints(dofelemid[oElemType.dof_type],
                                                                 nodes,
                                                                 nidxmaster,
                                                                 nidxslave,
                                                                 maxdist,
                                                                 eConstrType,
                                                                 uConstraintFlags);
               break;
          case ElementTypes::ELEM_SFNCON6H:
               pElem = SurfToNodeConstr<ShapeTria6H>::BuildJoints(dofelemid[oElemType.dof_type],
                                                                  nodes,
                                                                  nidxmaster,
                                                                  nidxslave,
                                                                  maxdist,
                                                                  eConstrType,
                                                                  uConstraintFlags);
               break;               
          case ElementTypes::ELEM_SFNCON8:
               pElem = SurfToNodeConstr<ShapeQuad8>::BuildJoints(dofelemid[oElemType.dof_type],
                                                                 nodes,
                                                                 nidxmaster,
                                                                 nidxslave,
                                                                 maxdist,
                                                                 eConstrType,
                                                                 uConstraintFlags);
               break;
          default:
               FEM_ASSERT(false);
          }

          if ((uConstraintFlags & CF_ELEM_DOF_PRE_ALLOCATED) && dofelemid[oElemType.dof_type] > edof[oElemType.dof_type].rows()) {
               throw std::runtime_error("dof_map.edof is not consistent with elements");
          }

          FEM_ASSERT(pElem != nullptr);

          rgElemBlocks.emplace_back(std::move(pElem));
     }
}
#endif

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
// PKG_ADD: autoload("FEM_MAT_ACCEL_LOAD", "__mboct_fem_pkg__.oct");
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
// PKG_ADD: autoload("FEM_DO_THERMAL", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_DO_STRUCTURAL", "__mboct_fem_pkg__.oct");

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
// PKG_DEL: autoload("FEM_MAT_ACCEL_LOAD", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_VEC_LOAD_CONSISTENT", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_VEC_LOAD_LUMPED", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_VEC_STRESS_CAUCH", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_VEC_STRAIN_TOTAL", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_SCA_STRESS_VMIS", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_CT_FIXED", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_CT_SLIDING", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_MAT_THERMAL_COND", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_MAT_HEAT_CAPACITY", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_DO_THERMAL", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_DO_STRUCTURAL", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_VEC_LOAD_THERMAL", "__mboct_fem_pkg__.oct", "remove");

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

     if (args.length() != 2) {
          print_usage();
          return retval;
     }

     const octave_scalar_map m_mesh = args(0).scalar_map_value();

     if (error_state) {
          return retval;
     }

     const auto it_nodes = m_mesh.seek("nodes");

     if (it_nodes == m_mesh.end()) {
          error("field mesh.nodes not found");
          return retval;
     }

     const Matrix nodes = m_mesh.contents(it_nodes).matrix_value();

     if (error_state) {
          return retval;
     }

     if (nodes.ndims() != 2 || nodes.columns() != 6) {
          error("mesh.nodes must be a Nx6 matrix");
          return retval;
     }

     const octave_scalar_map m_load_case = args(1).scalar_map_value();

     if (error_state) {
          return retval;
     }

     const auto it_locked_dof = m_load_case.seek("locked_dof");

     if (it_locked_dof == m_load_case.end()) {
          error("field load_case.locked_dof not found");
          return retval;
     }

     const boolNDArray locked_dof = m_load_case.contents(it_locked_dof).bool_array_value();

     if (error_state) {
          return retval;
     }

     const auto it_domain = m_load_case.seek("domain");

     DofMap::DomainType eDomain = DofMap::DO_STRUCTURAL;
     
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
         locked_dof.rows() != nodes.rows()) {
          error("size of load_case.locked_dof is not valid");
          return retval;
     }

     const auto iter_elements = m_mesh.seek("elements");

     if (iter_elements == m_mesh.end()) {
          error("field mesh.elements not found");
          return retval;
     }

     const octave_scalar_map m_elements = m_mesh.contents(iter_elements).scalar_map_value();

     if (error_state) {
          return retval;
     }

     boolNDArray dof_in_use(dim_vector(nodes.rows(), nodes.columns()), false);

     for (octave_idx_type i = 0; i < ElementTypes::iGetNumTypes(); ++i) {
          const auto& oElemType = ElementTypes::GetType(i);

          const auto iter_elem_type = m_elements.seek(oElemType.name);

          if (iter_elem_type == m_elements.end()) {
               continue;
          }

          switch (oElemType.type) {
          case ElementTypes::ELEM_ISO8:
          case ElementTypes::ELEM_ISO20:
          case ElementTypes::ELEM_PENTA15:
          case ElementTypes::ELEM_TET10H:
          case ElementTypes::ELEM_TET10: {
               const int32NDArray elnodes = m_elements.contents(iter_elem_type).int32_array_value();

               if (error_state) {
                    return retval;
               }

               for (octave_idx_type j = 0; j < elnodes.rows(); ++j) {
                    for (octave_idx_type k = 0; k < elnodes.columns(); ++k) {
                         const octave_idx_type idxnode = elnodes.xelem(j, k).value();

                         if (idxnode < 1 || idxnode > nodes.rows()) {
                              error("node index %Ld of element mesh.elements.%s(%Ld, %Ld) out of range %Ld:%Ld",
                                    static_cast<long long>(idxnode),
                                    oElemType.name,
                                    static_cast<long long>(j),
                                    static_cast<long long>(k),
                                    1LL,
                                    static_cast<long long>(nodes.rows()));
                              return retval;
                         }

                         for (octave_idx_type l = 0; l < 3; ++l) {
                              dof_in_use.xelem(idxnode - 1, l) = true;
                         }
                    }
               }
          } break;
          case ElementTypes::ELEM_BEAM2: {
               const octave_map m_beam2 = m_elements.contents(iter_elem_type).map_value();

               if (error_state) {
                    return retval;
               }

               const auto iter_nodes = m_beam2.seek("nodes");

               if (iter_nodes == m_beam2.end()) {
                    error("missing field mesh.elements.beam2.nodes");
                    return retval;
               }

               const Cell& ov_nodes = m_beam2.contents(iter_nodes);

               if (error_state) {
                    return retval;
               }

               for (octave_idx_type j = 0; j < ov_nodes.numel(); ++j) {
                    const int32NDArray elnodes = ov_nodes(j).int32_array_value();

                    if (error_state) {
                         return retval;
                    }

                    for (octave_idx_type l = 0; l < elnodes.numel(); ++l) {
                         const octave_idx_type idxnode = elnodes(l).value();

                         if (idxnode < 1 || idxnode > nodes.rows()) {
                              error("invalid node index for mesh.elements.%s(%Ld).nodes(%Ld)=%Ld",
                                    oElemType.name,
                                    static_cast<long long>(j),
                                    static_cast<long long>(l),
                                    static_cast<long long>(idxnode));
                              return retval;
                         }

                         for (octave_idx_type k = 0; k < 6; ++k) {
                              dof_in_use(idxnode - 1, k) = true;
                         }
                    }
               }
          } break;
          case ElementTypes::ELEM_RBE3: {
               const octave_map m_rbe3 = m_elements.contents(iter_elem_type).map_value();

               if (error_state) {
                    return retval;
               }

               const auto iter_nodesr = m_rbe3.seek("nodes");

               if (iter_nodesr == m_rbe3.end()) {
                    error("missing field mesh.elements.rbe3.nodes");
                    return retval;
               }

               const Cell& ov_nodesr = m_rbe3.contents(iter_nodesr);

               for (octave_idx_type j = 0; j < ov_nodesr.numel(); ++j) {
                    const int32NDArray elnodes = ov_nodesr(j).int32_array_value();

                    if (error_state) {
                         return retval;
                    }

                    if (!elnodes.numel()) {
                         error("invalid number of nodes for mesh.elements.%s(%Ld).nodes", oElemType.name, static_cast<long long>(j));
                         return retval;
                    }

                    const octave_idx_type idxnode = elnodes(0).value();

                    if (idxnode < 1 || idxnode > nodes.rows()) {
                         error("invalid node index for mesh.elements.%s(%Ld).nodes(1)=%Ld",
                               oElemType.name,
                               static_cast<long long>(j),
                               static_cast<long long>(idxnode));
                         return retval;
                    }

                    for (octave_idx_type k = 0; k < 6; ++k) {
                         dof_in_use(idxnode - 1, k) = true;
                    }
               }
          } break;

          default:
               continue;
          }
     }

     octave_idx_type icurrdof = 0;
     
     int32NDArray ndof(dim_vector(nodes.rows(), iNodeMaxDofIndex), -1);

     for (octave_idx_type i = 0; i < ndof.rows(); ++i) {
          for (octave_idx_type j = 0; j < ndof.columns(); ++j) {
               if (dof_in_use.xelem(i, j)) {
                    ndof.xelem(i, j) = locked_dof.xelem(i, j) ? 0 : ++icurrdof;
               }
          }
     }

     octave_map dof_map;

     dof_map.assign("ndof", octave_value(ndof));
     dof_map.assign("domain", octave_value(eDomain));

     int32NDArray idx_node(dim_vector(icurrdof, 1), -1);
     octave_idx_type icurrndof = 0;

     for (octave_idx_type i = 0; i < ndof.rows(); ++i) {
          for (octave_idx_type j = 0; j < ndof.columns(); ++j) {
               if (ndof.xelem(i, j).value() > 0) {
                    idx_node.xelem(icurrndof++) = ndof.xelem(i, j);
               }
          }
     }

     FEM_ASSERT(icurrndof == icurrdof);

     dof_map.assign("idx_node", octave_value(idx_node));

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
                    break;

               default:
                    continue;
               }

               const auto iter_etype = m_elements.seek(oElemType.name);

               if (iter_etype == m_elements.end()) {
                    continue;
               }

               const octave_map m_etype = m_elements.contents(iter_etype).map_value();

               if (error_state) {
                    return retval;
               }

               if (!m_etype.numel()) {
                    continue;
               }

               switch (oElemType.type) {
               case ElementTypes::ELEM_JOINT:
               case ElementTypes::ELEM_SFNCON4:
               case ElementTypes::ELEM_SFNCON6:
               case ElementTypes::ELEM_SFNCON6H:
               case ElementTypes::ELEM_SFNCON8: {
                    const char* const fn = oElemType.type == ElementTypes::ELEM_JOINT
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
                    case ElementTypes::ELEM_SFNCON8: {
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
                                   icurrconstr = ov_fn(j).rows();
                                   edof[CS_JOINT].elem_count++;
                                   break;
                              case ElementTypes::ELEM_SFNCON4:
                              case ElementTypes::ELEM_SFNCON6:
                              case ElementTypes::ELEM_SFNCON6H:
                              case ElementTypes::ELEM_SFNCON8:
#if HAVE_NLOPT == 1
                                   icurrconstr = SurfToNodeConstrBase::iGetNumDof(ov_constr, j);
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
                                   icurrconstr = ov_fn(j).rows();
                                   icurrjoints = 1;
                                   break;
                              case ElementTypes::ELEM_SFNCON4:
                              case ElementTypes::ELEM_SFNCON6:
                              case ElementTypes::ELEM_SFNCON6H:
                              case ElementTypes::ELEM_SFNCON8:
#if HAVE_NLOPT == 1
                                   icurrconstr = SurfToNodeConstrBase::iGetNumDof(ov_constr, j);
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
                                   edof[CS_RBE3].dof.xelem(j, l) = ++icurrdof;
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
          for (octave_idx_type i = 0; i < edof[k].dof.rows(); ++i) {
               for (octave_idx_type j = 0; j < edof[k].dof.columns(); ++j) {
                    const octave_idx_type idxedof = edof[k].dof.xelem(i, j).value();

                    if (idxedof > 0) {
                         idx_lambda(icurrlambda++) = idxedof;
                    }
               }
          }
     }

     idx_lambda.sort();

     FEM_ASSERT(icurrlambda == inumlambda);

     if (m_edof.nfields()) {
          dof_map.assign("edof", octave_value(m_edof));
     }

     if (inumlambda) {
          dof_map.assign("idx_lambda", octave_value(idx_lambda));
     }

     dof_map.assign("totdof", octave_value(icurrdof));

     retval.append(dof_map);

     return retval;
}

DEFUN_DLD(fem_pre_mesh_constr_surf_to_node, args, nargout,
          "-*- texinfo -*-\n"
          "@deftypefn {} @var{joints} = fem_pre_mesh_constr_surf_to_node(@var{nodes}, @var{elements})\n"
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

     if (nargin != 2) {
          print_usage();
          return retval;
     }

     const Matrix nodes(args(0).matrix_value());

     if (error_state) {
          return retval;
     }

     const octave_scalar_map elements(args(1).scalar_map_value());

     if (error_state) {
          return retval;
     }

     try {
          if (nodes.columns() != 6) {
               throw std::runtime_error("invalid number of columns for matrix nodes");
          }

          array<int32NDArray, DofMap::ELEM_TYPE_COUNT> edof;
          array<octave_idx_type, DofMap::ELEM_TYPE_COUNT> dofelemid = {0};

          vector<std::unique_ptr<ElementBlockBase> > rgElemBlocks;

          rgElemBlocks.reserve(2);

          for (octave_idx_type k = ElementTypes::ELEM_SFNCON4; k <= ElementTypes::ELEM_SFNCON8; ++k) {
               constexpr unsigned uFlags = SurfToNodeConstrBase::CF_IGNORE_NODES_OUT_OF_RANGE;
               const ElementTypes::TypeInfo& oElemType = ElementTypes::GetType(k);

               FEM_ASSERT(oElemType.type == k);

               SurfToNodeConstrBase::BuildJoints(nodes, elements, edof, dofelemid, oElemType, uFlags, rgElemBlocks);
          }

          octave_idx_type iNumElem = 0;

          for (const auto& oElemBlk: rgElemBlocks) {
               iNumElem += oElemBlk->iGetNumElem();
          }

          constexpr size_t nFields = 2;
          static const char* const rgFields[nFields] = {"C", "nodes"};
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

          if (error_state) {
               throw std::runtime_error("argument mesh must be a scalar struct");
          }

          const auto it_nodes = mesh.seek("nodes");

          if (it_nodes == mesh.end()) {
               throw std::runtime_error("missing field mesh.nodes in argument mesh");
          }

          const Matrix nodes(mesh.contents(it_nodes).matrix_value());

          if (error_state) {
               throw std::runtime_error("mesh.nodes must be a real matrix in argument mesh");
          }

          const auto it_elements = mesh.seek("elements");

          if (it_elements == mesh.end()) {
               throw std::runtime_error("missing field mesh.elements in argument mesh");
          }

          const octave_scalar_map elements(mesh.contents(it_elements).scalar_map_value());

          if (error_state) {
               throw std::runtime_error("mesh.elements must be a scalar struct in argument mesh");
          }

          const auto it_materials = mesh.seek("materials");

          if (it_materials == mesh.end()) {
               throw std::runtime_error("missing field mesh.materials in argument mesh");
          }

          const octave_scalar_map materials(mesh.contents(it_materials).scalar_map_value());

          if (error_state) {
               throw std::runtime_error("mesh.materials must be a scalar struct in argument mesh");
          }

          const auto it_material_data = mesh.seek("material_data");

          if (it_material_data == mesh.end()) {
               throw std::runtime_error("missing field mesh.material_data in argument mesh");
          }

          const octave_map material_data(mesh.contents(it_material_data).map_value());

          if (error_state) {
               throw std::runtime_error("mesh.material_data must be a struct array in argument mesh");
          }

          const octave_scalar_map dof_map(args(1).scalar_map_value());

          if (error_state) {
               throw std::runtime_error("argument dof_map must be a scalar struct");
          }

          const int32NDArray matrix_type(args(2).int32_array_value());

          if (error_state) {
               throw std::runtime_error("argument matrix_type must be an array of integers");
          }

          const octave_map load_case(nargin > 3 ? args(3).map_value() : octave_map());

          if (error_state) {
               throw std::runtime_error("argument load case must be a struct array");
          }

          const StrainField oRefStrain(load_case, nodes);

          const octave_scalar_map sol(nargin > 4 ? args(4).scalar_map_value() : octave_scalar_map());

          if (error_state) {
               throw std::runtime_error("argument sol must be scalar struct");
          }

          if (nodes.columns() != 6) {
               throw std::runtime_error("invalid number of columns for matrix mesh.nodes in argument mesh");
          }

          const auto iter_ndof = dof_map.seek("ndof");

          if (iter_ndof == dof_map.end()) {
               throw std::runtime_error("field \"ndof\" not found in argument dof_map");
          }

          const auto iter_domain = dof_map.seek("domain");

          if (iter_domain == dof_map.end()) {
               throw std::runtime_error("field \"domain\" not found in argument dof_map");
          }

          const auto eDomain = static_cast<DofMap::DomainType>(dof_map.contents(iter_domain).int_value());

          const int32NDArray ndof(dof_map.contents(iter_ndof).int32_array_value());

          if (error_state) {
               throw std::runtime_error("field dof_map.ndof must be an integer array in argument dof_map");
          }

          const auto iter_totdof = dof_map.seek("totdof");

          if (iter_totdof == dof_map.end()) {
               throw std::runtime_error("field \"totdof\" not found in argument dof_map");
          }

          const octave_idx_type inumdof = dof_map.contents(iter_totdof).int32_scalar_value();

          if (error_state) {
               throw std::runtime_error("field dof_map.totdof must be an scalar integer in argument dof_map");
          }

          if (ndof.rows() != nodes.rows() || ndof.columns() != DofMap::iGetNodeMaxDofIndex(eDomain)) {
               throw std::runtime_error("shape of dof_map.ndof is not valid");
          }

          for (octave_idx_type j = 0; j < ndof.columns(); ++j) {
               for (octave_idx_type i = 0; i < ndof.rows(); ++i) {
                    octave_idx_type idof = ndof(i, j).value();
                    if (idof > inumdof) {
                         throw std::runtime_error("invalid index in matrix dof_map.ndof in argument dof_map");
                    }
               }
          }

          const auto iter_U = sol.seek("def");

          const NDArray sol_U = (iter_U != sol.end()) ? sol.contents(iter_U).array_value() : NDArray();

          if (error_state) {
               throw std::runtime_error("field sol.def must be an array in argument sol");
          }

          const auto iter_edof = dof_map.seek("edof");

          array<int32NDArray, DofMap::ELEM_TYPE_COUNT> edof;
          array<octave_idx_type, DofMap::ELEM_TYPE_COUNT> dofelemid = {0};

          if (iter_edof != dof_map.end()) {
               const octave_scalar_map s_edof(dof_map.contents(iter_edof).scalar_map_value());

               if (error_state) {
                    throw std::runtime_error("dof_map.edof must be a scalar struct in argument dof_map");
               }

               static const struct DofEntries {
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

                         if (error_state) {
                              std::ostringstream os;
                              os << "dof_map.edof." << k->name << " must be an integer array in argument dof_map";
                              throw std::runtime_error(os.str());
                         }

                         if (a_edof.columns() < k->col_min || (k->col_max > 0 && a_edof.columns() > k->col_max)) {
                              std::ostringstream os;
                              os << "dof_map.edof." << k->name << " numer of columns = " << a_edof.columns() << " not in range [" << k->col_min << ":" << k->col_max << "] in argument dof_map";
                              throw std::runtime_error(os.str());
                         }

                         edof[k->type] = a_edof;
                    }
               }
          }

          for (auto k = std::begin(edof); k != std::end(edof); ++k) {
               for (octave_idx_type j = 0; j < k->columns(); ++j) {
                    for (octave_idx_type i = 0; i < k->rows(); ++i) {
                         if ((*k)(i, j).value() > inumdof) {
                              throw std::runtime_error("dof_map.edof dof index out of range in argument dof_map");
                         }
                    }
               }
          }

          const auto iterE = material_data.seek("E");
          const auto iternu = material_data.seek("nu");
          const auto iterC = material_data.seek("C");
          const auto iterRho = material_data.seek("rho");
          const auto iterAlpha = material_data.seek("alpha");
          const auto iterBeta = material_data.seek("beta");
          const auto iterGamma = material_data.seek("gamma");
          const auto iterk = material_data.seek("k");
          const auto itercp = material_data.seek("cp");
          
          if (iterC == material_data.end() && (iterE == material_data.end() || iternu == material_data.end())) {
               throw std::runtime_error("field \"C\" not found in mesh.material_data in argument mesh");
          }

          if (iterRho == material_data.end()) {
               throw std::runtime_error("field \"rho\" not found in mesh.material_data in argument mesh");
          }

          const Cell cellC = iterC != material_data.end() ? material_data.contents(iterC) : Cell();
          const Cell cellE = iterE != material_data.end() ? material_data.contents(iterE) : Cell();
          const Cell cellnu = iternu != material_data.end() ? material_data.contents(iternu) : Cell();
          const Cell cellRho = material_data.contents(iterRho);
          const Cell cellAlpha = iterAlpha != material_data.end() ? material_data.contents(iterAlpha) : Cell();
          const Cell cellBeta = iterBeta != material_data.end() ? material_data.contents(iterBeta) : Cell();
          const Cell cellGamma = iterGamma != material_data.end() ? material_data.contents(iterGamma) : Cell();
          const Cell cellk = iterk != material_data.end() ? material_data.contents(iterk) : Cell();
          const Cell cellcp = itercp != material_data.end() ? material_data.contents(itercp) : Cell();
          
          vector<Material> rgMaterials;

          rgMaterials.reserve(material_data.numel());

          Matrix C;

          for (octave_idx_type i = 0; i < material_data.numel(); ++i) {
               const bool buseC = iterC != material_data.end() && !cellC(i).isempty();
               double E, nu;

               if (buseC) {
                    C = cellC(i).matrix_value();

                    if (error_state || C.rows() != 6 || C.columns() != 6) {
                         throw std::runtime_error("size of constitutive matrix mesh.material_data.C is not valid in argument mesh");
                    }

                    if ((iterE != material_data.end() && !cellE(i).isempty()) || (iternu != material_data.end() && !cellnu(i).isempty())) {
                         throw std::runtime_error("redundant material properties in field material_data in argument mesh");
                    }
               } else {
                    if (iterE == material_data.end() || cellE(i).isempty() || iternu == material_data.end() || cellnu(i).isempty()) {
                         throw std::runtime_error("fields \"E\" and \"nu\" not found in mesh.material_data in argument mesh");
                    }

                    E = cellE(i).scalar_value();

                    if (error_state) {
                         throw std::runtime_error("field \"E\" in mesh.material_data must be a real scalar");
                    }

                    nu = cellnu(i).scalar_value();

                    if (error_state) {
                         throw std::runtime_error("field \"nu\" in mesh.material_data must be a real scalar");
                    }

                    if (iterC != material_data.end() && !cellC(i).isempty()) {
                         throw std::runtime_error("redundant material properties in field material_data in argument mesh");
                    }
               }

               const double rho = cellRho(i).scalar_value();

               if (error_state) {
                    throw std::runtime_error("mesh.material_data.rho is not a valid scalar in argument mesh");
               }

               const double alpha = iterAlpha != material_data.end() && !cellAlpha(i).isempty() ? cellAlpha(i).scalar_value() : 0.;

               if (error_state) {
                    throw std::runtime_error("mesh.material_data.alpha is not a valid scalar in argument mesh");
               }

               const double beta = iterBeta != material_data.end() && !cellBeta(i).isempty() ? cellBeta(i).scalar_value() : 0.;

               if (error_state) {
                    throw std::runtime_error("mesh.material_data.beta is not a valid scalar in argument mesh");
               }

               const double gamma = iterGamma != material_data.end() && !cellGamma(i).isempty() ? cellGamma(i).scalar_value() : 0.;

               const Matrix k = iterk != material_data.end() && !cellk(i).isempty() ? cellk(i).matrix_value() : Matrix(3, 3, 0.);

               if (k.rows() != 3 || k.columns() != 3) {
                    throw std::runtime_error("mesh.material_data.k must be a real 3x3 matrix");
               }

               const double cp = itercp != material_data.end() && !cellcp(i).isempty() ? cellcp(i).scalar_value() : 0.;

               if (buseC) {
                    rgMaterials.emplace_back(C, rho, alpha, beta, gamma, k, cp);
               } else {
                    rgMaterials.emplace_back(E, nu, rho, alpha, beta, gamma, k, cp);
               }
          }

          DofMap oDof(eDomain, ndof, edof, inumdof);

          array<bool, ElementTypes::iGetNumTypes()> rgElemUse;

          std::fill(std::begin(rgElemUse), std::end(rgElemUse), false);

          for (octave_idx_type i = 0; i < matrix_type.numel(); ++i) {
               const auto eMatType = static_cast<Element::FemMatrixType>(matrix_type(i).value());

               if ((eMatType & eDomain) != eDomain) {
                    throw std::runtime_error("matrix type is not valid for selected domain");
               }
               
               switch (eMatType) {
               case Element::MAT_STIFFNESS:
               case Element::MAT_STIFFNESS_SYM:
               case Element::MAT_STIFFNESS_SYM_L:
                    rgElemUse[ElementTypes::ELEM_RBE3] = true;
                    rgElemUse[ElementTypes::ELEM_JOINT] = true;
                    rgElemUse[ElementTypes::ELEM_SFNCON4] = true;
                    rgElemUse[ElementTypes::ELEM_SFNCON6] = true;
                    rgElemUse[ElementTypes::ELEM_SFNCON6H] = true;
                    rgElemUse[ElementTypes::ELEM_SFNCON8] = true;
                    // fall through
               case Element::MAT_MASS:
               case Element::MAT_MASS_SYM:
               case Element::MAT_MASS_SYM_L:
               case Element::MAT_MASS_LUMPED:
               case Element::MAT_DAMPING:
               case Element::MAT_DAMPING_SYM:
               case Element::MAT_DAMPING_SYM_L:
               case Element::SCA_TOT_MASS:
               case Element::VEC_INERTIA_M1:
               case Element::MAT_INERTIA_J:
               case Element::MAT_INERTIA_INV3:
               case Element::MAT_INERTIA_INV4:
               case Element::MAT_INERTIA_INV5:
               case Element::MAT_INERTIA_INV8:
               case Element::MAT_INERTIA_INV9:
               case Element::MAT_ACCEL_LOAD:
                    rgElemUse[ElementTypes::ELEM_BEAM2] = true;
                    // fall through
               case Element::VEC_STRESS_CAUCH:
               case Element::VEC_STRAIN_TOTAL:
               case Element::SCA_STRESS_VMIS:
               case Element::MAT_THERMAL_COND:
                    rgElemUse[ElementTypes::ELEM_THERM_CONV_ISO4] = true;
                    rgElemUse[ElementTypes::ELEM_THERM_CONV_QUAD8] = true;
                    rgElemUse[ElementTypes::ELEM_THERM_CONV_TRIA6] = true;
                    rgElemUse[ElementTypes::ELEM_THERM_CONV_TRIA6H] = true;
               case Element::MAT_HEAT_CAPACITY:
                    rgElemUse[ElementTypes::ELEM_ISO8] = true;
                    rgElemUse[ElementTypes::ELEM_ISO20] = true;
                    rgElemUse[ElementTypes::ELEM_PENTA15] = true;
                    rgElemUse[ElementTypes::ELEM_TET10H] = true;
                    rgElemUse[ElementTypes::ELEM_TET10] = true;
                    break;

               case Element::VEC_LOAD_CONSISTENT:
               case Element::VEC_LOAD_LUMPED:
                    if (load_case.numel() == 0) {
                         throw std::runtime_error("missing argument load_case for matrix_type == FEM_VEC_LOAD_*");
                    }

                    rgElemUse[ElementTypes::ELEM_PRESSURE_ISO4] = true;
                    rgElemUse[ElementTypes::ELEM_PRESSURE_TRIA6] = true;
                    rgElemUse[ElementTypes::ELEM_PRESSURE_TRIA6H] = true;
                    rgElemUse[ElementTypes::ELEM_PRESSURE_QUAD8] = true;
                    rgElemUse[ElementTypes::ELEM_STRUCT_FORCE] = true;
                    rgElemUse[ElementTypes::ELEM_JOINT] = true;
                    // Needed for thermal stress only
                    rgElemUse[ElementTypes::ELEM_ISO8] = true;
                    rgElemUse[ElementTypes::ELEM_ISO20] = true;
                    rgElemUse[ElementTypes::ELEM_PENTA15] = true;
                    rgElemUse[ElementTypes::ELEM_TET10H] = true;
                    rgElemUse[ElementTypes::ELEM_TET10] = true;                    
                    break;

               case Element::VEC_LOAD_THERMAL:
                    rgElemUse[ElementTypes::ELEM_THERM_CONV_ISO4] = true;
                    rgElemUse[ElementTypes::ELEM_THERM_CONV_QUAD8] = true;
                    rgElemUse[ElementTypes::ELEM_THERM_CONV_TRIA6] = true;
                    rgElemUse[ElementTypes::ELEM_THERM_CONV_TRIA6H] = true;
                    break;
                    
               default:
                    throw std::runtime_error("invalid value for argument matrix_type");
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
               case ElementTypes::ELEM_PENTA15:
               case ElementTypes::ELEM_TET10H:
               case ElementTypes::ELEM_TET10: {
                    const auto iter_elem = elements.seek(oElemType.name);

                    if (iter_elem == elements.end()) {
                         continue;
                    }

                    const int32NDArray elem_nodes = elements.contents(iter_elem).int32_array_value();

                    if (error_state) {
                         throw std::runtime_error(std::string("mesh.elements.") + oElemType.name
                                                  + " must be an array of integers in argument mesh");
                    }

                    if (elem_nodes.columns() < oElemType.min_nodes) {
                         throw std::runtime_error(std::string("invalid number of nodes for element type ") + oElemType.name
                                                  +  " in argument mesh");
                    }

                    if (oElemType.max_nodes > 0 && elem_nodes.columns() > oElemType.max_nodes) {
                         throw std::runtime_error(std::string("invalid number of nodes for element type ") + oElemType.name
                                                  + " in argument mesh");
                    }

                    for (octave_idx_type j = 0; j < elem_nodes.columns(); ++j) {
                         for (octave_idx_type i = 0; i < elem_nodes.rows(); ++i) {
                              octave_idx_type inode = elem_nodes(i, j);
                              if (inode < 1 || inode > nodes.rows()) {
                                   throw std::runtime_error(std::string("invalid node index for element type ")
                                                            + oElemType.name + " in argument mesh");
                              }
                         }
                    }

                    const auto iter_mat = materials.seek(oElemType.name);

                    if (iter_mat == materials.end()) {
                         throw std::runtime_error("material not defined for all element types in argument mesh");
                    }

                    const int32NDArray elem_mat = materials.contents(iter_mat).int32_array_value();

                    if (elem_mat.rows() != elem_nodes.rows()) {
                         throw std::runtime_error("invalid number of rows for matrix mesh.materials in argument mesh");
                    }

                    for (octave_idx_type i = 0; i < elem_mat.rows(); ++i) {
                         octave_idx_type imaterial = elem_mat(i);
                         if (imaterial <= 0 || imaterial > material_data.numel()) {
                              throw std::runtime_error("invalid index in matrix mesh.materials in argument mesh");
                         }
                    }

                    switch (oElemType.type) {
                    case ElementTypes::ELEM_ISO8:
                         rgElemBlocks.emplace_back(new ElementBlock<Iso8>(oElemType.type, elem_nodes, nodes, 3, elem_mat, rgMaterials, oRefStrain));
                         break;

                    case ElementTypes::ELEM_ISO20:
                         rgElemBlocks.emplace_back(new ElementBlock<Iso20>(oElemType.type, elem_nodes, nodes, 3, elem_mat, rgMaterials, oRefStrain));
                         break;

                    case ElementTypes::ELEM_PENTA15:
                         rgElemBlocks.emplace_back(new ElementBlock<Penta15>(oElemType.type, elem_nodes, nodes, 3, elem_mat, rgMaterials, oRefStrain));
                         break;

                    case ElementTypes::ELEM_TET10H:
                         rgElemBlocks.emplace_back(new ElementBlock<Tet10h>(oElemType.type, elem_nodes, nodes, 3, elem_mat, rgMaterials, oRefStrain));
                         break;

                    case ElementTypes::ELEM_TET10:
                         rgElemBlocks.emplace_back(new ElementBlock<Tet10>(oElemType.type, elem_nodes, nodes, 3, elem_mat, rgMaterials, oRefStrain));
                         break;

                    default:
                         throw std::runtime_error("invalid element type");
                    }
               } break;
               case ElementTypes::ELEM_BEAM2:
               case ElementTypes::ELEM_RBE3:
               case ElementTypes::ELEM_JOINT: {
                    const auto iter_elem = elements.seek(oElemType.name);

                    if (iter_elem == elements.end()) {
                         continue;
                    }

                    const octave_map s_elem(elements.contents(iter_elem).map_value());

                    if (error_state) {
                         throw std::runtime_error("mesh.elements."s + oElemType.name + " must be an struct array in argument mesh");
                    }

                    const auto iter_nodes = s_elem.seek("nodes");

                    if (iter_nodes == s_elem.end()) {
                         throw std::runtime_error("missing field mesh.elements."s + oElemType.name + ".nodes in argument mesh");
                    }

                    const Cell ov_nodes(s_elem.contents(iter_nodes));

                    FEM_ASSERT(ov_nodes.numel() == s_elem.numel());

                    Cell ov_material, ov_section, ov_e2;

                    if (oElemType.type == ElementTypes::ELEM_BEAM2) {
                         const auto iter_mat = s_elem.seek("material");

                         if (iter_mat == s_elem.end()) {
                              throw std::runtime_error("missing field mesh.elements.beam2.material in argument mesh");
                         }

                         ov_material = s_elem.contents(iter_mat);

                         FEM_ASSERT(ov_material.numel() == s_elem.numel());

                         const auto iter_sect = s_elem.seek("section");

                         if (iter_sect == s_elem.end()) {
                              throw std::runtime_error("missing field mesh.elements.beam2.section in argument mesh");
                         }

                         ov_section = s_elem.contents(iter_sect);

                         FEM_ASSERT(ov_section.numel() == s_elem.numel());

                         const auto iter_e2 = s_elem.seek("e2");

                         if (iter_e2 == s_elem.end()) {
                              throw std::runtime_error("missing field mesh.elements.beam2.e2 in argument mesh");
                         }

                         ov_e2 = s_elem.contents(iter_e2);

                         FEM_ASSERT(ov_e2.numel() == s_elem.numel());
                    }

                    Cell ov_C;

                    if (oElemType.type == ElementTypes::ELEM_JOINT) {
                         const auto iter_C = s_elem.seek("C");

                         if (iter_C == s_elem.end()) {
                              throw std::runtime_error("missing field mesh.elements.joints.C in argument mesh");
                         }

                         ov_C = s_elem.contents(iter_C);

                         FEM_ASSERT(ov_C.numel() == s_elem.numel());
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
                    case ElementTypes::ELEM_RBE3:
                         pElem.reset(new ElementBlock<ElemRBE3>(oElemType.type, s_elem.numel()));
                         break;
                    case ElementTypes::ELEM_JOINT:
                         pElem.reset(new ElementBlock<ElemJoint>(oElemType.type, s_elem.numel()));
                         break;
                    default:
                         FEM_ASSERT(false);
                    }

                    for (octave_idx_type i = 0; i < s_elem.numel(); ++i) {
                         const int32NDArray elem_nodes(ov_nodes(i).int32_array_value());

                         if (error_state) {
                              throw std::runtime_error("mesh.elements"s + oElemType.name + ".nodes must be an array of integers in argument mesh");
                         }

                         if (elem_nodes.columns() < oElemType.min_nodes) {
                              throw std::runtime_error("invalid number of nodes for element type "s + oElemType.name + " in argument mesh");
                         }

                         if (oElemType.max_nodes > 0 && elem_nodes.columns() > oElemType.max_nodes) {
                              throw std::runtime_error("invalid number of nodes for element type "s + oElemType.name + " in argument mesh");
                         }

                         if (elem_nodes.rows() != 1) {
                              throw std::runtime_error("invalid number of rows in node matrix for element type "s + oElemType.name + " in argument mesh");
                         }

                         for (octave_idx_type j = 0; j < elem_nodes.columns(); ++j) {
                              octave_idx_type inode = elem_nodes(j);

                              if (inode < 1 || inode > nodes.rows()) {
                                   throw std::runtime_error("invalid node index for element type "s + oElemType.name + " in argument mesh");
                              }
                         }

                         Matrix X(6, elem_nodes.columns());

                         for (octave_idx_type j = 0; j < elem_nodes.columns(); ++j) {
                              for (octave_idx_type k = 0; k < X.rows(); ++k) {
                                   X.xelem(k, j) = nodes.xelem(elem_nodes(j).value() - 1, k);
                              }
                         }

                         switch (oElemType.type) {
                         case ElementTypes::ELEM_BEAM2: {
                              const octave_idx_type imaterial(ov_material(i).int_value() - 1);

                              if (error_state) {
                                   throw std::runtime_error("mesh.elements.beam2.material must be a scalar integer");
                              }

                              if (imaterial < 0 || static_cast<size_t>(imaterial) >= rgMaterials.size()) {
                                   throw std::runtime_error("mesh.elements.beam2.material out of range");
                              }

                              const Material* const material = &rgMaterials[imaterial];

                              const octave_scalar_map m_section(ov_section(i).scalar_map_value());

                              if (error_state) {
                                   throw std::runtime_error("mesh.elements.beam2.section must be a scalar struct");
                              }

                              auto iterA = m_section.seek("A");

                              if (iterA == m_section.end()) {
                                   throw std::runtime_error("missing field mesh.elements.beam2.section.A");
                              }

                              auto iterAy = m_section.seek("Ay");

                              if (iterAy == m_section.end()) {
                                   throw std::runtime_error("missing field mesh.elements.beam2.section.Ay");
                              }

                              auto iterAz = m_section.seek("Az");

                              if (iterAz == m_section.end()) {
                                   throw std::runtime_error("missing field mesh.elements.beam2.section.Az");
                              }

                              auto iterIt = m_section.seek("It");

                              if (iterIt == m_section.end()) {
                                   throw std::runtime_error("missing field mesh.elements.beam2.section.It");
                              }

                              auto iterIy = m_section.seek("Iy");

                              if (iterIy == m_section.end()) {
                                   throw std::runtime_error("missing field mesh.elements.beam2.section.Iy");
                              }

                              auto iterIz = m_section.seek("Iz");

                              if (iterIz == m_section.end()) {
                                   throw std::runtime_error("missing field mesh.elements.beam2.section.Iz");
                              }

                              BeamCrossSection section;

                              section.A = m_section.contents(iterA).scalar_value();

                              if (error_state) {
                                   throw std::runtime_error("mesh.elements.beam2.section.A must be a real scalar");
                              }

                              section.Ay = m_section.contents(iterAy).scalar_value();

                              if (error_state) {
                                   throw std::runtime_error("mesh.elements.beam2.section.Ay must be a real scalar");
                              }
                              section.Az = m_section.contents(iterAz).scalar_value();

                              if (error_state) {
                                   throw std::runtime_error("mesh.elements.beam2.section.Az must be a real scalar");
                              }

                              section.It = m_section.contents(iterIt).scalar_value();

                              if (error_state) {
                                   throw std::runtime_error("mesh.elements.beam2.section.It must be a real scalar");
                              }

                              section.Iy = m_section.contents(iterIy).scalar_value();

                              if (error_state) {
                                   throw std::runtime_error("mesh.elements.beam2.section.Iy must be a real scalar");
                              }

                              section.Iz = m_section.contents(iterIz).scalar_value();

                              if (error_state) {
                                   throw std::runtime_error("mesh.elements.beam2.section.Iz must be a real scalar");
                              }

                              const ColumnVector e2(ov_e2(i).column_vector_value());

                              if (error_state) {
                                   throw std::runtime_error("mesh.elements.beam2.e2 must be a column vector");
                              }

                              if (e2.rows() != 3) {
                                   throw std::runtime_error("mesh.elements.beam2.e2 must be 3x1 matrix");
                              }

                              pElem->Insert<ElemBeam2>(i + 1, X, material, elem_nodes, section, e2);
                         } break;
                         case ElementTypes::ELEM_RBE3: {
                              RowVector weight;

                              if (ov_weight.numel() > i) {
                                   weight = ov_weight(i).row_vector_value();

                                   if (error_state) {
                                        throw std::runtime_error("mesh.elements.rbe3.weight must be a row vector in argument mesh");
                                   }
                              } else {
                                   weight.resize(elem_nodes.numel(), 1.);
                              }

                              if (weight.numel() != elem_nodes.numel() - 1) {
                                   throw std::runtime_error("numel(mesh.elements.rbe3.weight) does not match numel(mesh.elements.rbe3.nodes) - 1 in argument mesh");
                              }

                              pElem->Insert<ElemRBE3>(++dofelemid[oElemType.dof_type], X, nullptr, elem_nodes, weight);
                         } break;
                         case ElementTypes::ELEM_JOINT: {
                              const Matrix C(ov_C(i).matrix_value());

                              if (error_state) {
                                   throw std::runtime_error("mesh.elements.joints.C must be a real matrix in argument mesh");
                              }

                              if (C.rows() < 1 || C.rows() > edof[oElemType.dof_type].columns() || C.columns() != 6 * elem_nodes.columns() || C.rows() > C.columns()) {
                                   throw std::runtime_error("invalid size for field elements.joints.C");
                              }

                              Matrix U(C.rows(), load_case.numel(), 0.); // By default displacement is set to zero

                              const auto iter_joints = load_case.seek("joints");

                              if (iter_joints != load_case.end()) {
                                   const Cell ov_joints = load_case.contents(iter_joints);

                                   for (octave_idx_type k = 0; k < load_case.numel(); ++k) {
                                        if (ov_joints(k).isempty()) {
                                             continue;
                                        }

                                        const octave_map s_joints(ov_joints(k).map_value());

                                        if (error_state) {
                                             throw std::runtime_error("load_case.joints must be an struct array in argument load_case");
                                        }

                                        const auto iter_U = s_joints.seek("U");

                                        if (iter_U == s_elem.end()) {
                                             throw std::runtime_error("missing field load_case.joints.U in argument load case");
                                        }

                                        const Cell ov_U(s_joints.contents(iter_U));

                                        if (ov_U.numel() != s_elem.numel()) {
                                             throw std::runtime_error("load_case.joints must have the same size like mesh.elements.joints in argument load case");
                                        }

                                        const ColumnVector Uk(ov_U(i).column_vector_value());

                                        if (error_state || Uk.rows() != C.rows()) {
                                             throw std::runtime_error("load_case.joints.U must be a real column vector of the same number of rows like mesh.elements.joints.C in argument load_case");
                                        }

                                        for (octave_idx_type l = 0; l < C.rows(); ++l) {
                                             U.xelem(l, k) = Uk.xelem(l);
                                        }
                                   }
                              }

                              pElem->Insert<ElemJoint>(++dofelemid[oElemType.dof_type], X, nullptr, elem_nodes, C, U);
                         } break;
                         default:
                              FEM_ASSERT(false);
                         }
                    }

                    if (oElemType.dof_type != DofMap::ELEM_NODOF && dofelemid[oElemType.dof_type] > edof[oElemType.dof_type].rows()) {
                         throw std::runtime_error("dof_map.edof is not consistent with elements in argument dof_map");
                    }

                    rgElemBlocks.emplace_back(std::move(pElem));
               } break;
               case ElementTypes::ELEM_SFNCON4:
               case ElementTypes::ELEM_SFNCON6:
               case ElementTypes::ELEM_SFNCON6H:
               case ElementTypes::ELEM_SFNCON8: {
#if HAVE_NLOPT == 1
                    constexpr unsigned uFlags = SurfToNodeConstrBase::CF_ELEM_DOF_PRE_ALLOCATED;
                    SurfToNodeConstrBase::BuildJoints(nodes, elements, edof, dofelemid, oElemType, uFlags, rgElemBlocks);
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
               case ElementTypes::ELEM_PRESSURE_TRIA6:
                    InsertPressureElem<PressureLoadTria6>(oElemType.type, nodes, load_case, oElemType.name, oElemType.max_nodes, rgElemBlocks);
                    break;
               case ElementTypes::ELEM_PRESSURE_TRIA6H:
                    InsertPressureElem<PressureLoadTria6H>(oElemType.type, nodes, load_case, oElemType.name, oElemType.max_nodes, rgElemBlocks);
                    break;                    
               case ElementTypes::ELEM_STRUCT_FORCE: {
                    const auto iter_loads = load_case.seek("loads");
                    const auto iter_loaded_nodes = load_case.seek("loaded_nodes");

                    if (iter_loads != load_case.end() && iter_loaded_nodes != load_case.end()) {
                         const Cell cell_loads = load_case.contents(iter_loads);
                         const Cell cell_loaded_nodes = load_case.contents(iter_loaded_nodes);

                         FEM_ASSERT(cell_loads.numel() == cell_loaded_nodes.numel());
                         FEM_ASSERT(cell_loads.numel() == load_case.numel());

                         octave_idx_type iNumForces = 0;

                         for (octave_idx_type i = 0; i < cell_loads.numel(); ++i) {
                              if (cell_loads(i).numel()) {
                                   ++iNumForces;
                              }
                         }

                         if (iNumForces) {
                              std::unique_ptr<ElementBlock<StructForce> > pElem(new ElementBlock<StructForce>(oElemType.type, iNumForces));

                              for (octave_idx_type i = 0; i < cell_loads.numel(); ++i) {
                                   if (cell_loads(i).numel()) {
                                        const Matrix loads = cell_loads(i).matrix_value();

                                        if (error_state) {
                                             throw std::runtime_error("field load_case.loads must be a real matrix in argument load_case");
                                        }

                                        if (loads.columns() != 3 && loads.columns() != 6) {
                                             throw std::runtime_error("load_case.loads must be a n x 3 or n x 6 matrix in argument load_case");
                                        }

                                        const int32NDArray loaded_nodes = cell_loaded_nodes(i).int32_array_value();

                                        if (error_state) {
                                             throw std::runtime_error("field load_case.loaded_nodes must be an integer matrix in argument load_case");
                                        }

                                        if (loaded_nodes.columns() != 1 || loaded_nodes.rows() != loads.rows()) {
                                             throw std::runtime_error("load_case.loaded_nodes must be a column vector with the same number of rows like load_case.loads");
                                        }

                                        Matrix X(3, loaded_nodes.rows());

                                        for (octave_idx_type l = 0; l < X.columns(); ++l) {
                                             for (octave_idx_type m = 0; m < X.rows(); ++m) {
                                                  octave_idx_type inode = loaded_nodes(l).value() - 1;

                                                  if (inode < 0 || inode >= nodes.rows()) {
                                                       throw std::runtime_error("node index out of range in load_case.pressure.elements in argument load_case");
                                                  }

                                                  X.xelem(m, l) = nodes.xelem(inode, m);
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
                    InsertThermalConvElem<ThermalConvectionIso4>(oElemType.type, nodes, elements, load_case, oElemType.name, oElemType.max_nodes, rgElemBlocks);
                    break;
               case ElementTypes::ELEM_THERM_CONV_QUAD8:
                    InsertThermalConvElem<ThermalConvectionQuad8>(oElemType.type, nodes, elements, load_case, oElemType.name, oElemType.max_nodes, rgElemBlocks);
                    break;
               case ElementTypes::ELEM_THERM_CONV_TRIA6:
                    InsertThermalConvElem<ThermalConvectionTria6>(oElemType.type, nodes, elements, load_case, oElemType.name, oElemType.max_nodes, rgElemBlocks);
                    break;
               case ElementTypes::ELEM_THERM_CONV_TRIA6H:
                    InsertThermalConvElem<ThermalConvectionTria6H>(oElemType.type, nodes, elements, load_case, oElemType.name, oElemType.max_nodes, rgElemBlocks);
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
                    throw std::runtime_error("dof_map.edof is not consistent with mesh.elements in argument dof_map");
               }
          }

          octave_idx_type iMaxWorkSpaceSize = 0;

          for (octave_idx_type j = 0; j < matrix_type.numel(); ++j) {
               const Element::FemMatrixType eMatType = static_cast<Element::FemMatrixType>(matrix_type(j).value());
               octave_idx_type iWorkSpaceSize = 0;

               for (auto i = rgElemBlocks.cbegin(); i != rgElemBlocks.cend(); ++i) {
                    iWorkSpaceSize += (*i)->iGetWorkSpaceSize(eMatType);
               }

               iMaxWorkSpaceSize = std::max(iMaxWorkSpaceSize, iWorkSpaceSize);
          }

          MatrixAss oMatAss(iMaxWorkSpaceSize);
          MeshInfo oMeshInfo;
          MatrixAss::MatrixInfo oMatInfo;
          bool bMatInfo = false;

          for (octave_idx_type i = 0; i < matrix_type.numel(); ++i) {
               const Element::FemMatrixType eMatType = static_cast<Element::FemMatrixType>(matrix_type.xelem(i).value());

               switch (eMatType) {
               case Element::MAT_STIFFNESS:
               case Element::MAT_STIFFNESS_SYM:
               case Element::MAT_STIFFNESS_SYM_L:
               case Element::MAT_MASS:
               case Element::MAT_MASS_SYM:
               case Element::MAT_MASS_SYM_L:
               case Element::MAT_MASS_LUMPED:
               case Element::MAT_DAMPING:
               case Element::MAT_DAMPING_SYM:
               case Element::MAT_DAMPING_SYM_L:
               case Element::VEC_LOAD_CONSISTENT:
               case Element::VEC_LOAD_LUMPED:
               case Element::MAT_ACCEL_LOAD:
               case Element::MAT_THERMAL_COND:
               case Element::MAT_HEAT_CAPACITY:
               case Element::VEC_LOAD_THERMAL: {
                    oMatAss.Reset(eMatType, oMatInfo);

                    for (auto j = rgElemBlocks.cbegin(); j != rgElemBlocks.cend(); ++j) {
                         const bool bNeedMatInfo = (*j)->bNeedMatrixInfo(eMatType);

                         if (!bMatInfo && bNeedMatInfo) {
                              oMatAss.UpdateMatrixInfo();
                              oMatInfo = oMatAss.GetMatrixInfo();
                              bMatInfo = true;
                         }

                         FEM_TRACE("i=" << i << " beta=" << oMatInfo.beta << "\nalpha=" << oMatInfo.alpha << "\n");

                         (*j)->Assemble(oMatAss, oMeshInfo, oDof, eMatType);
                    }

                    oMatAss.Finish();
                    retval.append(oMatAss.Assemble(oDof, load_case.numel()));
               } break;
               case Element::SCA_TOT_MASS: {
                    double dMass = 0.;

                    for (auto j = rgElemBlocks.cbegin(); j != rgElemBlocks.cend(); ++j) {
                         dMass += (*j)->dGetMass();
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
                    const octave_idx_type iNumModes = sol_U.ndims() == 3 ? sol_U.dim3() : 0;
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
                              throw std::runtime_error("argument sol is not optional for selected matrix type in argument matrix_type");
                         }

                         if (sol_U.rows() != nodes.rows() || sol_U.columns() != nodes.columns()) {
                              throw std::runtime_error("dimensions of argument sol.def does not match dimension of mesh.nodes in argument mesh");
                         }

                         if (sol_U.ndims() != 3) {
                              throw std::runtime_error("sol.def must be a three dimensional array in argument sol");
                         }
                    }

                    NDArray mat(mat_dim, 0.);

                    for (auto j = rgElemBlocks.cbegin(); j != rgElemBlocks.cend(); ++j) {
                         (*j)->PostProcElem(mat, eMatType, sol_U);
                    }

                    switch (eMatType) {
                    case Element::MAT_INERTIA_J:
                         FEM_ASSERT(mat.ndims() == 2);
                         FEM_ASSERT(mat.rows() == 3);
                         FEM_ASSERT(mat.columns() == 3);

                         for (octave_idx_type i = 1; i < mat.rows(); ++i) {
                              for (octave_idx_type j = 0; j < i; ++j) {
                                   mat.xelem(i, j) = mat.xelem(j, i);
                              }
                         }
                         break;

                    default:
                         break;
                    }

                    retval.append(mat);
               } break;
               case Element::VEC_STRESS_CAUCH:
               case Element::VEC_STRAIN_TOTAL:
               case Element::SCA_STRESS_VMIS: {
                    if (nargin <= 4) {
                         throw std::runtime_error("argument sol is not optional for selected matrix type in argument matrix_type");
                    }

                    if (sol_U.rows() != nodes.rows() || sol_U.columns() != nodes.columns()) {
                         throw std::runtime_error("dimensions of argument sol.def does not match dimension of mesh.nodes in argument mesh");
                    }

                    if (sol_U.ndims() < 2) {
                         throw std::runtime_error("sol.def must be a three dimensional array in argument sol");
                    }

                    constexpr octave_idx_type iNumStress = 6;
                    const octave_idx_type iNumLoads = sol_U.ndims() >= 3 ? sol_U.dim3() : 1;

                    octave_scalar_map s_StressStrain, s_StressStrainm, s_vmis;

                    for (octave_idx_type j = 0; j < ElementTypes::iGetNumTypes(); ++j) {
                         const ElementTypes::TypeInfo& oElemType = ElementTypes::GetType(j);

                         if (!rgElemUse[oElemType.type]) {
                              continue;
                         }

                         switch (oElemType.type) {
                         case ElementTypes::ELEM_ISO8:
                         case ElementTypes::ELEM_ISO20:
                         case ElementTypes::ELEM_PENTA15:
                         case ElementTypes::ELEM_TET10H:
                         case ElementTypes::ELEM_TET10: {
                              const auto iter_elem = elements.seek(oElemType.name);

                              if (iter_elem == elements.end()) {
                                   continue;
                              }

                              const int32NDArray elem_nodes = elements.contents(iter_elem).int32_array_value();

                              const Element::FemMatrixType eMatTypeStressStrain = eMatType == Element::SCA_STRESS_VMIS ? Element::VEC_STRESS_CAUCH : eMatType;

                              NDArray oElemStressStrain(dim_vector(elem_nodes.rows(), elem_nodes.columns(), iNumStress, iNumLoads), 0.);

                              for (auto k = rgElemBlocks.cbegin(); k != rgElemBlocks.cend(); ++k) {
                                   if ((*k)->GetElementType() == oElemType.type) {
                                        (*k)->PostProcElem(oElemStressStrain, eMatTypeStressStrain, sol_U);
                                   }
                              }

                              intNDArray<octave_idx_type> iStressStrain(dim_vector(nodes.rows(), 1), 0);

                              for (octave_idx_type k = 0; k < oElemStressStrain.dim1(); ++k) {
                                   for (octave_idx_type l = 0; l < oElemStressStrain.dim2(); ++l) {
                                        const octave_idx_type inode = elem_nodes.xelem(k, l).value() - 1;
                                        ++iStressStrain.xelem(inode);
                                   }
                              }

                              NDArray oNodalStressStrain(dim_vector(nodes.rows(), iNumStress, iNumLoads), 0.);

                              for (octave_idx_type n = 0; n < iNumLoads; ++n) {
                                   for (octave_idx_type m = 0; m < iNumStress; ++m) {
                                        for (octave_idx_type k = 0; k < oElemStressStrain.dim1(); ++k) {
                                             for (octave_idx_type l = 0; l < oElemStressStrain.dim2(); ++l) {
                                                  const octave_idx_type inode = elem_nodes.xelem(k, l).value() - 1;
                                                  oNodalStressStrain.xelem(inode, m, n) += oElemStressStrain.xelem(k, l, m + n * iNumStress) / iStressStrain.xelem(inode);
                                             }
                                        }
                                   }
                              }

                              NDArray oContStressStrain(oElemStressStrain.dims());

                              for (octave_idx_type n = 0; n < iNumLoads; ++n) {
                                   for (octave_idx_type m = 0; m < iNumStress; ++m) {
                                        for (octave_idx_type k = 0; k < oElemStressStrain.dim1(); ++k) {
                                             for (octave_idx_type l = 0; l < oElemStressStrain.dim2(); ++l) {
                                                  const octave_idx_type inode = elem_nodes.xelem(k, l).value() - 1;
                                                  oContStressStrain.xelem(k, l, m + n * iNumStress) = oNodalStressStrain.xelem(inode, m, n);
                                             }
                                        }
                                   }
                              }

                              if (eMatType == Element::SCA_STRESS_VMIS) {
                                   NDArray vmis(dim_vector(elem_nodes.rows(), elem_nodes.columns(), iNumLoads));

                                   for (octave_idx_type n = 0; n < iNumLoads; ++n) {
                                        for (octave_idx_type l = 0; l < oElemStressStrain.dim2(); ++l) {
                                             for (octave_idx_type k = 0; k < oElemStressStrain.dim1(); ++k) {
                                                  const octave_idx_type ioffset = n * iNumStress;

                                                  const double tauxx = oContStressStrain.xelem(k, l, ioffset);
                                                  const double tauyy = oContStressStrain.xelem(k, l, ioffset + 1);
                                                  const double tauzz = oContStressStrain.xelem(k, l, ioffset + 2);
                                                  const double tauxy = oContStressStrain.xelem(k, l, ioffset + 3);
                                                  const double tauyz = oContStressStrain.xelem(k, l, ioffset + 4);
                                                  const double tauzx = oContStressStrain.xelem(k, l, ioffset + 5);

                                                  vmis.xelem(k, l, n) = sqrt(tauxx * tauxx + tauyy * tauyy + tauzz * tauzz
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
               default:
                    throw std::runtime_error("invalid value for argument matrix_type");
               }
          }

          octave_scalar_map mat_info;

          mat_info.assign("alpha", oMatInfo.alpha);
          mat_info.assign("beta", oMatInfo.beta);

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

DEFINE_GLOBAL_CONSTANT(Element, MAT_ACCEL_LOAD, "acceleration load matrix (e.g. for gravity loads)")
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
DEFINE_GLOBAL_CONSTANT(Element, MAT_STIFFNESS_SYM, "upper triangular part of the stiffness matrix")
DEFINE_GLOBAL_CONSTANT(Element, MAT_STIFFNESS_SYM_L, "lower triangular part of the stiffness matrix")
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
DEFINE_GLOBAL_CONSTANT(DofMap, DO_STRUCTURAL, "structural domain")
DEFINE_GLOBAL_CONSTANT(DofMap, DO_THERMAL, "thermal domain")
