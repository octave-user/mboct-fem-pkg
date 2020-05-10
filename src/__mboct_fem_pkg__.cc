// Copyright (C) 2018(-2020) Reinhard <octave-user@a1.net>
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

#if HAVE_NLOPT == 1
#include <nlopt.h>
#endif

#ifdef DEBUG
#define OCTAVE_ENABLE_BOUNDS_CHECK
#define FEM_ASSERT(expr) if (!(expr)) {                         \
                std::cerr << "assertion " << #expr              \
                          << " failed at " << __FILE__ << ":"   \
                          << __LINE__ << ":"                    \
                          << __FUNCTION__ << std::endl;         \
                asm("int3");                                    \
                assert(expr);                                   \
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
        enum ElementType {
                ELEM_RBE3 = 0,
                ELEM_JOINT,
                ELEM_TYPE_COUNT,
                ELEM_NODOF = -1
        };

        explicit DofMap(const int32NDArray& ndof, const array<int32NDArray, ELEM_TYPE_COUNT>& edof, octave_idx_type totdof)
                :ndof(ndof), edof(edof), totdof(totdof) {
        }

        explicit DofMap(const int32NDArray& ndof, octave_idx_type totdof)
                :ndof(ndof), totdof(totdof) {
        }

        octave_idx_type GetNodeDofIndex(octave_idx_type inode, octave_idx_type idof) const {
                return ndof(inode, idof);
        }

        octave_idx_type GetElemDofIndex(ElementType eElemType, octave_idx_type ielem, octave_idx_type idof) const {
                return edof[eElemType](ielem, idof).value();
        }

        octave_idx_type iGetNumDof() const {
                return totdof;
        }

private:
        int32NDArray ndof;
        array<int32NDArray, ELEM_TYPE_COUNT> edof;
        octave_idx_type totdof;
};

class MatrixAss;

class MeshInfo {
public:
        enum InfoType {
                JACOBIAN_DET = 0,
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

                static const char szInfoName[INFO_COUNT][5] = {
                        "detJ"
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
        Material(const Matrix& C, double rho, double alpha, double beta)
                :rho(rho), alpha(alpha), beta(beta), C(C) {
                FEM_ASSERT(C.rows() == 6);
                FEM_ASSERT(C.columns() == 6);
        }

        Material(double E, double nu, double rho, double alpha, double beta)
                :rho(rho), alpha(alpha), beta(beta), C(6, 6) {

                C(0,0) = (E*(1-nu))/((1-2*nu)*(nu+1));
                C(1,0) = (E*nu)/((1-2*nu)*(nu+1));
                C(2,0) = (E*nu)/((1-2*nu)*(nu+1));
                C(3,0) = 0;
                C(4,0) = 0;
                C(5,0) = 0;
                C(0,1) = (E*nu)/((1-2*nu)*(nu+1));
                C(1,1) = (E*(1-nu))/((1-2*nu)*(nu+1));
                C(2,1) = (E*nu)/((1-2*nu)*(nu+1));
                C(3,1) = 0;
                C(4,1) = 0;
                C(5,1) = 0;
                C(0,2) = (E*nu)/((1-2*nu)*(nu+1));
                C(1,2) = (E*nu)/((1-2*nu)*(nu+1));
                C(2,2) = (E*(1-nu))/((1-2*nu)*(nu+1));
                C(3,2) = 0;
                C(4,2) = 0;
                C(5,2) = 0;
                C(0,3) = 0;
                C(1,3) = 0;
                C(2,3) = 0;
                C(3,3) = (E/(nu+1))/2.0e+0;
                C(4,3) = 0;
                C(5,3) = 0;
                C(0,4) = 0;
                C(1,4) = 0;
                C(2,4) = 0;
                C(3,4) = 0;
                C(4,4) = (E/(nu+1))/2.0e+0;
                C(5,4) = 0;
                C(0,5) = 0;
                C(1,5) = 0;
                C(2,5) = 0;
                C(3,5) = 0;
                C(4,5) = 0;
                C(5,5) = (E/(nu+1))/2.0e+0;
        }

        const Matrix& LinearElasticity() const {
                return C;
        }

        double Density() const { return rho; }
        double AlphaDamping() const { return alpha; }
        double BetaDamping() const { return beta; }
private:
        double rho, alpha, beta;
        Matrix C;
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
                r(iEvalPnt, iDirection) = ri;
        }

        void SetWeight(octave_idx_type iEvalPnt, double alphai) {
                alpha(iEvalPnt) = alphai;
        }

        double dGetPosition(octave_idx_type iEvalPnt, octave_idx_type iDirection) const {
                return r(iEvalPnt, iDirection);
        }

        double dGetWeight(octave_idx_type iEvalPnt) const {
                return alpha(iEvalPnt);
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

        enum MatrixType: unsigned {
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
                SCA_STRESS_VMIS = 0x1000u,
        };

        Element(octave_idx_type id, const Matrix& X, const Material* material, const int32NDArray& nodes)
                :id(id), X(X), material(material), nodes(nodes) {

                FEM_ASSERT(X.columns() == nodes.numel());
        }

        virtual ~Element() {
        }

        virtual void Assemble(MatrixAss& mat, MeshInfo& info, const DofMap& dof, MatrixType eMatType) const=0;

        virtual void PostProcElem(NDArray& mat, MatrixType eMatType, const NDArray& U) const {
        }

        virtual octave_idx_type iGetWorkSpaceSize(MatrixType eMatType) const=0;

        virtual double dGetMass() const {
                return 0.;
        }

        virtual bool bNeedMatrixInfo(MatrixType eMatType) const {
                return false;
        }

        virtual void Extract(octave_idx_type& idx, octave_map& sElem) const {
        }

protected:
        octave_idx_type id;
        Matrix X;
        const Material* material;
        int32NDArray nodes;
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
                        if (ridx(i).value() == cidx(i).value()) {
                                diagA(ridx(i)) += data(i);
                        }
                }

                constexpr double INIT_MIN = std::numeric_limits<double>::max();
                constexpr double INIT_MAX = -std::numeric_limits<double>::max();

                double minA = INIT_MIN;
                double maxA = INIT_MAX;

                for (octave_idx_type i = 0; i < nnz; ++i) {
                        const double absA = fabs(diagA(i));

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
                                Insert(Ke(i, j), r(i), c(j));
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

        void Reset(Element::MatrixType eMatTypeCurr, const MatrixInfo& oMatInfo) {
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
                
                data(nnz) = d;
                ridx(nnz) = r;
                cidx(nnz) = c;
                ++nnz;
        }

        Element::MatrixType eMatType;
        octave_idx_type nnz;
        int32NDArray ridx, cidx;
        ColumnVector data;
        MatrixInfo info;
};

class ElemJoint: public Element
{
public:
        ElemJoint(octave_idx_type id, const Matrix& X, const Material* material, const int32NDArray& nodes, const Matrix& C, const Matrix& U)
                :Element(id, X, material, nodes), C(C), U(U) {
                FEM_ASSERT(C.columns() == nodes.numel() * 6);
                FEM_ASSERT(C.rows() <= C.columns());
                FEM_ASSERT(C.rows() >= 1);
                FEM_ASSERT(U.rows() == C.rows());
        }

        virtual void Assemble(MatrixAss& mat, MeshInfo& info, const DofMap& dof, const MatrixType eMatType) const {
                switch (eMatType) {
                case MAT_STIFFNESS:
                case MAT_STIFFNESS_SYM:
                case MAT_STIFFNESS_SYM_L: {
                        int32NDArray ndofidx(dim_vector(nodes.numel() * 6, 1), -1);

                        for (octave_idx_type inode = 0; inode < nodes.numel(); ++inode) {
                                for (octave_idx_type idof = 0; idof < 6; ++idof) {
                                        ndofidx(inode * 6 + idof) = dof.GetNodeDofIndex(nodes(inode).value() - 1, idof);
                                }
                        }

                        int32NDArray edofidx(dim_vector(C.rows(), 1));

                        for (octave_idx_type idof = 0; idof < edofidx.numel(); ++idof) {
                                edofidx(idof) = dof.GetElemDofIndex(DofMap::ELEM_JOINT, id - 1, idof);
                        }

                        const double beta = mat.GetMatrixInfo().beta;


                        for (octave_idx_type j = 0; j < C.columns(); ++j) {
                                for (octave_idx_type i = 0; i < C.rows(); ++i) {
                                        const double Cij = beta * C(i, j);
                                        mat.Insert(Cij, ndofidx(j), edofidx(i));
                                        mat.Insert(Cij, edofidx(i), ndofidx(j));
                                }
                        }
                } break;
                case VEC_LOAD_CONSISTENT:
                case VEC_LOAD_LUMPED: {
                        int32NDArray edofidx(dim_vector(C.rows(), 1));

                        for (octave_idx_type idof = 0; idof < edofidx.numel(); ++idof) {
                                edofidx(idof) = dof.GetElemDofIndex(DofMap::ELEM_JOINT, id - 1, idof);
                        }

                        const double beta = mat.GetMatrixInfo().beta;

                        for (octave_idx_type j = 0; j < U.columns(); ++j) {
                                for (octave_idx_type i = 0; i < U.rows(); ++i) {
                                        mat.Insert(beta * U(i, j), edofidx(i), j + 1);
                                }
                        }
                } break;
                default:
                        break;
                }
        }

        virtual octave_idx_type iGetWorkSpaceSize(MatrixType eMatType) const {
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

        virtual bool bNeedMatrixInfo(Element::MatrixType eMatType) const {
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
        Matrix C, U;
};

class ElemRBE3: public Element
{
public:
        ElemRBE3(octave_idx_type id, const Matrix& X, const Material* material, const int32NDArray& nodes, const RowVector& omega)
                :Element(id, X, material, nodes),
                 omega(omega) {

                FEM_ASSERT(X.rows() == 6);
                FEM_ASSERT(X.columns() > 1);
                FEM_ASSERT(omega.numel() == nodes.numel() - 1);
        }

        virtual void Assemble(MatrixAss& mat, MeshInfo& info, const DofMap& dof, const MatrixType eMatType) const {
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
                                ndofidx(inode * 6 + idof) = dof.GetNodeDofIndex(nodes(inode).value() - 1, idof);
                        }
                }

                int32NDArray edofidx(dim_vector(6, 1));

                for (octave_idx_type idof = 0; idof < edofidx.rows(); ++idof) {
                        edofidx(idof) = dof.GetElemDofIndex(DofMap::ELEM_RBE3, id - 1, idof);
                }

                Matrix xi(3, nodes.numel() - 1);

                for (octave_idx_type j = 1; j < nodes.numel(); ++j) {
                        for (octave_idx_type i = 0; i < 3; ++i) {
                                xi(i, j - 1) = X(i, j) - X(i, 0);
                        }
                }

                FEM_TRACE("xi=[\n" << xi << "];\n");

                Matrix S(xi.columns() * 6, 6);

                for (octave_idx_type k = 0; k < xi.columns(); ++k) {
                        for (octave_idx_type j = 0; j < 6; ++j) {
                                for (octave_idx_type i = 0; i < 6; ++i) {
                                        const bool alpha = j < 3 || ndofidx((k + 1) * 6 + j).value() >= 0;
                                        S(6 * k + i, j) = alpha * (i == j);
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

                        S(6 * k + 0, 4) =  xi(2, k);
                        S(6 * k + 0, 5) = -xi(1, k);
                        S(6 * k + 1, 3) = -xi(2, k);
                        S(6 * k + 1, 5) =  xi(0, k);
                        S(6 * k + 2, 3) =  xi(1, k);
                        S(6 * k + 2, 4) = -xi(0, k);
                }

                FEM_TRACE("S=[\n" << S << "];\n");

                double Lc2 = 0.;

                for (octave_idx_type k = 0; k < xi.columns(); ++k) {
                        double norm_xik = 0.;

                        for (octave_idx_type i = 0; i < 3; ++i) {
                                norm_xik += xi(i, k) * xi(i, k);
                        }

                        Lc2 += sqrt(norm_xik);
                }

                Lc2 /= xi.columns();
                Lc2 *= Lc2;

                FEM_TRACE("Lc2=" << Lc2 << ";\n");

                ColumnVector W(xi.columns() * 6);

                for (octave_idx_type k = 0; k < xi.columns(); ++k) {
                        const double omegak = omega(k);

                        for (octave_idx_type i = 0; i < 6; ++i) {
                                W(k * 6 + i) = omegak * (i < 3 ? 1. : Lc2);
                        }
                }

                Matrix STWS(6, 6, 0.);

                for (octave_idx_type j = 0; j < 6; ++j) {
                        for (octave_idx_type i = 0; i < 6; ++i) {
                                for (octave_idx_type k = 0; k < S.rows(); ++k) {
                                        STWS(i, j) += S(k, i) * W(k) * S(k, j);
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
                                                Bijl += W(l * 6 + i) * S(l * 6 + i, k) * X(k, j);
                                        }

                                        B(l * 6 + i, j) = Bijl;
                                }
                        }
                }

                FEM_TRACE("f=[\n" << (B - DiagMatrix(W) * S * X) << "];\n");
                FEM_TRACE("W=[\n" << W << "];\n");
                FEM_TRACE("B=[\n" << B << "];\n");

                const double beta = mat.GetMatrixInfo().beta;

                for (octave_idx_type i = 0; i < 6; ++i) {
                        mat.Insert(-beta, ndofidx(i), edofidx(i));
                        mat.Insert(-beta, edofidx(i), ndofidx(i));
                }

                for (octave_idx_type j = 0; j < 6; ++j) {
                        for (octave_idx_type i = 0; i < xi.columns() * 6; ++i) {
                                const double Bij = beta * B(i, j);
                                mat.Insert(Bij, ndofidx(i + 6), edofidx(j));
                                mat.Insert(Bij, edofidx(j), ndofidx(i + 6));
                        }
                }
        }

        virtual octave_idx_type iGetWorkSpaceSize(MatrixType eMatType) const {
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

        virtual bool bNeedMatrixInfo(Element::MatrixType eMatType) const {
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

class Element3D: public Element
{
public:
        Element3D(octave_idx_type id, const Matrix& X, const Material* material, const int32NDArray& nodes)
                :Element(id, X, material, nodes) {

                FEM_ASSERT(X.rows() == 3);
        }

        virtual void Assemble(MatrixAss& mat, MeshInfo& info, const DofMap& dof, MatrixType eMatType) const {
                void (Element3D::*pFunc)(Matrix&, MeshInfo&, MatrixType) const;

                const octave_idx_type iNumDof = iGetNumDof();

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

                case MAT_ACCEL_LOAD:
                        pFunc = &Element3D::AccelerationLoad;
                        iNumRows = iNumDof;
                        iNumCols = 3;
                        break;

                default:
                        return;
                }

                int32NDArray dofidx(dim_vector(iNumDof, 1), 0);

                for (octave_idx_type inode = 0; inode < nodes.numel(); ++inode) {
                        for (octave_idx_type idof = 0; idof < 3; ++idof) {
                                dofidx(inode * 3 + idof) = dof.GetNodeDofIndex(nodes(inode).value() - 1, idof);
                        }
                }

                Matrix Ae(iNumRows, iNumCols, 0.);

                (this->*pFunc)(Ae, info, eMatType);

                switch (eMatType) {
                case MAT_ACCEL_LOAD: {
                        int32NDArray dofidxcol(dim_vector(iNumCols, 1));

                        for (octave_idx_type i = 0; i < iNumCols; ++i) {
                                dofidxcol(i) = i + 1;
                        }

                        mat.Insert(Ae, dofidx, dofidxcol);
                } break;
                default:
                        mat.Insert(Ae, dofidx, dofidx);
                }
        }

        octave_idx_type iGetNumDof() const {
                return nodes.numel() * 3;
        }

        virtual octave_idx_type iGetWorkSpaceSize(MatrixType eMatType) const {
                switch (eMatType) {
                case MAT_MASS:
                case MAT_STIFFNESS:
                case MAT_DAMPING: {
                        const octave_idx_type iNumDof = iGetNumDof();

                        return iNumDof * iNumDof;
                }
                case MAT_MASS_SYM:
                case MAT_MASS_SYM_L:
                case MAT_STIFFNESS_SYM:
                case MAT_STIFFNESS_SYM_L:
                case MAT_DAMPING_SYM:
                case MAT_DAMPING_SYM_L: {
                        const octave_idx_type iNumDof = iGetNumDof();

                        return (iNumDof + 1) * iNumDof / 2;
                }
                case MAT_MASS_LUMPED:
                        return iGetNumDof();

                case MAT_ACCEL_LOAD:
                        return iGetNumDof() * 3;

                default:
                        return 0;
                }
        }

        virtual void PostProcElem(NDArray& mat, MatrixType eMatType, const NDArray& U) const {
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
                                rv(j) = oIntegRule.dGetPosition(i, j);
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
        virtual const IntegrationRule& GetIntegrationRule(MatrixType eMatType) const=0;
        virtual double Jacobian(const ColumnVector& rv, Matrix& J) const=0;
        virtual void StrainMatrix(const ColumnVector& rv, const Matrix& invJ, Matrix& B) const=0;
        virtual void DispInterpMatrix(const ColumnVector& rv, Matrix& H) const=0;
        virtual Matrix InterpGaussToNodal(MatrixType eMatType, const Matrix& taun) const=0;

        void AddMeshInfo(MeshInfo& info, const IntegrationRule& oIntegRule, double detJ) const {
                info.Add(MeshInfo::JACOBIAN_DET, detJ);
        }

        void StiffnessMatrix(Matrix& Ke, MeshInfo& info, MatrixType eMatType) const {
                const IntegrationRule& oIntegRule = GetIntegrationRule(eMatType);
                const octave_idx_type iNumDof = iGetNumDof();
                const octave_idx_type iNumDir = oIntegRule.iGetNumDirections();
                ColumnVector rv(iNumDir);
                const Matrix& C = material->LinearElasticity();
                const octave_idx_type iNumStrains = C.rows();

                FEM_ASSERT(C.rows() == C.columns());

                Matrix J(iNumDir, iNumDir), B(iNumStrains, iNumDof), CB(iNumStrains, iNumDof);

                for (octave_idx_type i = 0; i < oIntegRule.iGetNumEvalPoints(); ++i) {
                        const double alpha = oIntegRule.dGetWeight(i);

                        for (octave_idx_type j = 0; j < iNumDir; ++j) {
                                rv(j) = oIntegRule.dGetPosition(i, j);
                        }

                        const double detJ = Jacobian(rv, J);

                        AddMeshInfo(info, oIntegRule, detJ);

                        StrainMatrix(rv, J.inverse(), B);

                        for (octave_idx_type l = 0; l < iNumStrains; ++l) {
                                for (octave_idx_type m = 0; m < iNumDof; ++m) {
                                        double CBlm = 0.;

                                        for (octave_idx_type n = 0; n < iNumStrains; ++n) {
                                                CBlm += detJ * alpha * C(l, n) * B(n, m);
                                        }

                                        FEM_ASSERT(std::isfinite(CBlm));

                                        CB(l, m) = CBlm;
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
                                                Kelm += B(n, l) * CB(n, m);
                                        }

                                        FEM_ASSERT(std::isfinite(Kelm));

                                        Ke(l, m) += Kelm;
                                }
                        }
                }

                for (octave_idx_type i = 1; i < iNumDof; ++i) {
                        for (octave_idx_type j = 0; j < i; ++j) {
                                Ke(i, j) = Ke(j, i);
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

        void MassMatrix(Matrix& Me, MeshInfo& info, MatrixType eMatType) const {
                const IntegrationRule& oIntegRule = GetIntegrationRule(eMatType);
                const octave_idx_type iNumDof = iGetNumDof();
                const octave_idx_type iNumDir = oIntegRule.iGetNumDirections();
                ColumnVector rv(iNumDir);
                const double rho = material->Density();
                const octave_idx_type iNumDisp = X.rows();

                Matrix J(iNumDir, iNumDir), H(iNumDisp, iNumDof);

                for (octave_idx_type i = 0; i < oIntegRule.iGetNumEvalPoints(); ++i) {
                        const double alpha = oIntegRule.dGetWeight(i);

                        for (octave_idx_type j = 0; j < iNumDir; ++j) {
                                rv(j) = oIntegRule.dGetPosition(i, j);
                        }

                        const double detJ = Jacobian(rv, J);

                        AddMeshInfo(info, oIntegRule, detJ);

                        DispInterpMatrix(rv, H);

                        for (octave_idx_type l = 0; l < iNumDof; ++l) {
                                for (octave_idx_type m = l; m < iNumDof; ++m) {
                                        double Melm = 0.;

                                        for (octave_idx_type n = 0; n < iNumDisp; ++n) {
                                                Melm += H(n, l) * H(n, m);
                                        }

                                        FEM_ASSERT(std::isfinite(Melm));

                                        Me(l, m) += Melm * alpha * rho * detJ;
                                }
                        }
                }

                for (octave_idx_type i = 1; i < iNumDof; ++i) {
                        for (octave_idx_type j = 0; j < i; ++j) {
                                Me(i, j) = Me(j, i);
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

        void DampingMatrix(Matrix& De, MeshInfo& info, MatrixType eMatType) const {
                const double alpha = material->AlphaDamping();

                if (alpha) {
                        MassMatrix(De, info, static_cast<MatrixType>(MAT_MASS | (eMatType & MAT_SYM_MASK)));

                        for (octave_idx_type j = 0; j < De.columns(); ++j) {
                                for (octave_idx_type i = 0; i < De.rows(); ++i) {
                                        De(i, j) *= alpha;
                                }
                        }
                }

                const double beta = material->BetaDamping();

                if (beta) {
                        Matrix Ke(De.rows(), De.columns(), 0.);

                        StiffnessMatrix(Ke, info, static_cast<MatrixType>(MAT_STIFFNESS | (eMatType & MAT_SYM_MASK)));

                        for (octave_idx_type j = 0; j < De.columns(); ++j) {
                                for (octave_idx_type i = 0; i < De.rows(); ++i) {
                                        De(i, j) += beta * Ke(i, j);
                                }
                        }
                }
        }

        void AccelerationLoad(Matrix& C1, MeshInfo& info, MatrixType eMatType) const {
                const IntegrationRule& oIntegRule = GetIntegrationRule(eMatType);
                const octave_idx_type iNumDof = iGetNumDof();
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
                                rv(j) = oIntegRule.dGetPosition(i, j);
                        }

                        const double detJ = Jacobian(rv, J);

                        AddMeshInfo(info, oIntegRule, detJ);

                        DispInterpMatrix(rv, H);

                        const double dm = alpha * rho * detJ;

                        for (octave_idx_type k = 0; k < iNumDisp; ++k) {
                                for (octave_idx_type j = 0; j < iNumDof; ++j) {
                                        C1(j, k) += H(k, j) * dm;
                                }
                        }
                }
        }

        void InertiaMoment1(NDArray& S, MatrixType eMatType) const {
                const IntegrationRule& oIntegRule = GetIntegrationRule(eMatType);
                const octave_idx_type iNumDof = iGetNumDof();
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
                                rv(j) = oIntegRule.dGetPosition(i, j);
                        }

                        const double detJ = Jacobian(rv, J);;

                        DispInterpMatrix(rv, H);

                        for (octave_idx_type l = 0; l < S.rows(); ++l) {
                                double fil = 0.;

                                for (octave_idx_type n = 0; n < iNumNodes; ++n) {
                                        for (octave_idx_type m = 0; m < iNumDisp; ++m) {
                                                fil += H(l, iNumDisp * n + m) * X(m, n);
                                        }
                                }

                                S(l) += fil * alpha * rho * detJ;
                        }
                }
        }

        void InertiaMatrix(NDArray& Inv7, MatrixType eMatType) const {
                const IntegrationRule& oIntegRule = GetIntegrationRule(eMatType);
                const octave_idx_type iNumDof = iGetNumDof();
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
                                rv(j) = oIntegRule.dGetPosition(i, j);
                        }

                        const double detJ = Jacobian(rv, J);

                        DispInterpMatrix(rv, H);

                        for (octave_idx_type l = 0; l < iNumDisp; ++l) {
                                fi(l) = 0.;

                                for (octave_idx_type n = 0; n < iNumNodes; ++n) {
                                        for (octave_idx_type m = 0; m < iNumDisp; ++m) {
                                                fi(l) += H(l, iNumDisp * n + m) * X(m, n);
                                        }
                                }
                        }

                        const double dmi = alpha * rho * detJ;

                        Inv7(0, 0) += (fi(1) * fi(1) + fi(2) * fi(2)) * dmi;
                        Inv7(0, 1) -= (fi(0) * fi(1)) * dmi;
                        Inv7(0, 2) -= (fi(0) * fi(2)) * dmi;
                        Inv7(1, 1) += (fi(0) * fi(0) + fi(2) * fi(2)) * dmi;
                        Inv7(1, 2) -= (fi(1) * fi(2)) * dmi;
                        Inv7(2, 2) += (fi(0) * fi(0) + fi(1) * fi(1)) * dmi;
                }
        }

        void InertiaInv3(NDArray& Inv3, MatrixType eMatType, const NDArray& U) const {
                const IntegrationRule& oIntegRule = GetIntegrationRule(eMatType);
                const octave_idx_type iNumDof = iGetNumDof();
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
                                rv(j) = oIntegRule.dGetPosition(i, j);
                        }

                        const double detJ = Jacobian(rv, J);

                        DispInterpMatrix(rv, H);

                        const double dmi = alpha * rho * detJ;

                        for (octave_idx_type j = 0; j < U.dim3(); ++j) {
                                for (octave_idx_type l = 0; l < iNumDisp; ++l) {
                                        double Uil = 0.;

                                        for (octave_idx_type m = 0; m < iNumDisp; ++m) {
                                                for (octave_idx_type n = 0; n < iNumNodes; ++n) {
                                                        Uil += H(l, iNumDisp * n + m) * U(nodes(n).value() - 1, m, j);
                                                }
                                        }

                                        Inv3(l, j) += dmi * Uil;
                                }
                        }
                }
        }

        void InertiaInv4(NDArray& Inv4, MatrixType eMatType, const NDArray& U) const {
                const IntegrationRule& oIntegRule = GetIntegrationRule(eMatType);
                const octave_idx_type iNumDof = iGetNumDof();
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
                                rv(j) = oIntegRule.dGetPosition(i, j);
                        }

                        const double detJ = Jacobian(rv, J);

                        DispInterpMatrix(rv, H);

                        const double dmi = alpha * rho * detJ;

                        for (octave_idx_type j = 0; j < U.dim3(); ++j) {
                                for (octave_idx_type l = 0; l < iNumDisp; ++l) {
                                        Ui(l) = 0.;
                                        fi(l) = 0.;

                                        for (octave_idx_type m = 0; m < iNumDisp; ++m) {
                                                for (octave_idx_type n = 0; n < iNumNodes; ++n) {
                                                        const double Hlmn = H(l, iNumDisp * n + m);
                                                        Ui(l) += Hlmn * U(nodes(n).value() - 1, m, j);
                                                        fi(l) += Hlmn * X(m, n);
                                                }
                                        }
                                }

                                Inv4(0, j) += (fi(1) * Ui(2) - Ui(1) * fi(2)) * dmi;
                                Inv4(1, j) += (Ui(0) * fi(2) - fi(0) * Ui(2)) * dmi;
                                Inv4(2, j) += (fi(0) * Ui(1) - Ui(0) * fi(1)) * dmi;
                        }
                }
        }

        void InertiaInv5(NDArray& Inv5, MatrixType eMatType, const NDArray& U) const {
                const IntegrationRule& oIntegRule = GetIntegrationRule(eMatType);
                const octave_idx_type iNumDof = iGetNumDof();
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
                                rv(j) = oIntegRule.dGetPosition(i, j);
                        }

                        const double detJ = Jacobian(rv, J);

                        DispInterpMatrix(rv, H);

                        dmi(i) = alpha * rho * detJ;

                        for (octave_idx_type j = 0; j < U.dim3(); ++j) {
                                for (octave_idx_type l = 0; l < iNumDisp; ++l) {
                                        for (octave_idx_type m = 0; m < iNumDisp; ++m) {
                                                for (octave_idx_type n = 0; n < iNumNodes; ++n) {
                                                        Ui(l, j, i) += H(l, iNumDisp * n + m) * U(nodes(n).value() - 1, m, j);
                                                }
                                        }
                                }
                        }
                }

                for (octave_idx_type i = 0; i < iNumGauss; ++i) {
                        for (octave_idx_type j = 0; j < U.dim3(); ++j) {
                                for (octave_idx_type k = 0; k < U.dim3(); ++k) {
                                        Inv5(0, k, j) += dmi(i) * (Ui(1, j, i) * Ui(2, k, i) - Ui(2, j, i) * Ui(1, k, i));
                                        Inv5(1, k, j) += dmi(i) * (Ui(2, j, i) * Ui(0, k, i) - Ui(0, j, i) * Ui(2, k, i));
                                        Inv5(2, k, j) += dmi(i) * (Ui(0, j, i) * Ui(1, k, i) - Ui(1, j, i) * Ui(0, k, i));
                                }
                        }
                }
        }

        void InertiaInv8(NDArray& Inv8, MatrixType eMatType, const NDArray& U) const {
                const IntegrationRule& oIntegRule = GetIntegrationRule(eMatType);
                const octave_idx_type iNumDof = iGetNumDof();
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
                                rv(j) = oIntegRule.dGetPosition(i, j);
                        }

                        const double detJ = Jacobian(rv, J);

                        DispInterpMatrix(rv, H);

                        const double dmi = alpha * rho * detJ;

                        for (octave_idx_type j = 0; j < U.dim3(); ++j) {
                                for (octave_idx_type l = 0; l < iNumDisp; ++l) {
                                        Ui(l) = 0.;
                                        fi(l) = 0.;

                                        for (octave_idx_type m = 0; m < iNumDisp; ++m) {
                                                for (octave_idx_type n = 0; n < iNumNodes; ++n) {
                                                        const double Hlmn = H(l, iNumDisp * n + m);
                                                        Ui(l) += Hlmn * U(nodes(n).value() - 1, m, j);
                                                        fi(l) += Hlmn * X(m, n);
                                                }
                                        }
                                }

                                Inv8(0,0,j) += (fi(2)*Ui(2)+fi(1)*Ui(1))*dmi;
                                Inv8(1,0,j) += -fi(0)*Ui(1)*dmi;
                                Inv8(2,0,j) += -fi(0)*Ui(2)*dmi;
                                Inv8(0,1,j) += -Ui(0)*fi(1)*dmi;
                                Inv8(1,1,j) += (fi(2)*Ui(2)+fi(0)*Ui(0))*dmi;
                                Inv8(2,1,j) += -fi(1)*Ui(2)*dmi;
                                Inv8(0,2,j) += -Ui(0)*fi(2)*dmi;
                                Inv8(1,2,j) += -Ui(1)*fi(2)*dmi;
                                Inv8(2,2,j) += (fi(1)*Ui(1)+fi(0)*Ui(0))*dmi;
                        }
                }
        }

        void InertiaInv9(NDArray& Inv9, MatrixType eMatType, const NDArray& U) const {
                const IntegrationRule& oIntegRule = GetIntegrationRule(eMatType);
                const octave_idx_type iNumDof = iGetNumDof();
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
                                rv(j) = oIntegRule.dGetPosition(i, j);
                        }

                        const double detJ = Jacobian(rv, J);

                        DispInterpMatrix(rv, H);

                        dmi(i) = alpha * rho * detJ;

                        for (octave_idx_type j = 0; j < U.dim3(); ++j) {
                                for (octave_idx_type l = 0; l < iNumDisp; ++l) {
                                        for (octave_idx_type m = 0; m < iNumDisp; ++m) {
                                                for (octave_idx_type n = 0; n < iNumNodes; ++n) {
                                                        Ui(l, j, i) += H(l, iNumDisp * n + m) * U(nodes(n).value() - 1, m, j);
                                                }
                                        }
                                }
                        }
                }

                for (octave_idx_type i = 0; i < iNumGauss; ++i) {
                        for (octave_idx_type j = 0; j < U.dim3(); ++j) {
                                for (octave_idx_type k = 0; k < U.dim3(); ++k) {
                                        const octave_idx_type jk = j + k * Inv9.dim3();

                                        Inv9(0, 0, jk) += (-Ui(2,j,i)*Ui(2,k,i)-Ui(1,j,i)*Ui(1,k,i))*dmi(i);
                                        Inv9(1, 0, jk) += Ui(0,j,i)*Ui(1,k,i)*dmi(i);
                                        Inv9(2, 0, jk) += Ui(0,j,i)*Ui(2,k,i)*dmi(i);
                                        Inv9(0, 1, jk) += Ui(0,k,i)*Ui(1,j,i)*dmi(i);
                                        Inv9(1, 1, jk) += (-Ui(2,j,i)*Ui(2,k,i)-Ui(0,j,i)*Ui(0,k,i))*dmi(i);
                                        Inv9(2, 1, jk) += Ui(1,j,i)*Ui(2,k,i)*dmi(i);
                                        Inv9(0, 2, jk) += Ui(0,k,i)*Ui(2,j,i)*dmi(i);
                                        Inv9(1, 2, jk) += Ui(1,k,i)*Ui(2,j,i)*dmi(i);
                                        Inv9(2, 2, jk) += (-Ui(1,j,i)*Ui(1,k,i)-Ui(0,j,i)*Ui(0,k,i))*dmi(i);
                                }
                        }
                }
        }

        void StressNodalElem(NDArray& taun, MatrixType eMatType, const NDArray& U) const {
                const IntegrationRule& oIntegRule = GetIntegrationRule(eMatType);
                const octave_idx_type iNumGauss = oIntegRule.iGetNumEvalPoints();
                const octave_idx_type iNumDof = iGetNumDof();
                const octave_idx_type iNumDir = oIntegRule.iGetNumDirections();
                const octave_idx_type iNumLoads = U.ndims() >= 3 ? U.dim3() : 1;
                const octave_idx_type iNumNodes = nodes.numel();

                ColumnVector rv(iNumDir);
                const Matrix& C = material->LinearElasticity();
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

                Matrix J(iNumDir, iNumDir), B(iNumStrains, iNumDof), CB(iNumStrains, iNumDof);
                ColumnVector Ue(iNumDof);
                Matrix taug(iNumGauss, iNumStrains * iNumLoads);

                for (octave_idx_type i = 0; i < iNumGauss; ++i) {
                        for (octave_idx_type j = 0; j < iNumDir; ++j) {
                                rv(j) = oIntegRule.dGetPosition(i, j);
                        }

                        Jacobian(rv, J);

                        StrainMatrix(rv, J.inverse(), B);

                        for (octave_idx_type j = 0; j < iNumStrains; ++j) {
                                for (octave_idx_type k = 0; k < iNumDof; ++k) {
                                        double CBjk = 0;

                                        for (octave_idx_type l = 0; l < iNumStrains; ++l) {
                                                CBjk += C(j, l) * B(l, k);
                                        }

                                        CB(j, k) = CBjk;
                                }
                        }

                        for (octave_idx_type l = 0; l < iNumLoads; ++l) {
                                for (octave_idx_type j = 0; j < iNumNodes; ++j) {
                                        for (octave_idx_type k = 0; k < 3; ++k) {
                                                FEM_ASSERT(nodes(j).value() > 0);
                                                FEM_ASSERT(nodes(j).value() <= U.dim1());
                                                Ue(3 * j + k) = U(nodes(j).value() - 1, k, l);
                                        }
                                }

                                for (octave_idx_type j = 0; j < iNumStrains; ++j) {
                                        double taugj = 0;

                                        for (octave_idx_type k = 0; k < iNumDof; ++k) {
                                                taugj += CB(j, k) * Ue(k);
                                        }

                                        taug(i, l * iNumStrains + j) = taugj;
                                }
                        }
                }

                const Matrix tauen = InterpGaussToNodal(eMatType, taug);

                FEM_ASSERT(tauen.rows() == iNumNodes);
                FEM_ASSERT(tauen.columns() == taug.columns());

                for (octave_idx_type k = 0; k < iNumLoads; ++k) {
                        for (octave_idx_type j = 0; j < iNumStrains; ++j) {
                                for (octave_idx_type i = 0; i < iNumNodes; ++i) {
                                        taun(id - 1, i, j + k * iNumStrains) = tauen(i, k * iNumStrains + j);
                                }
                        }
                }
        }
};

class Iso8: public Element3D
{
public:
        Iso8(octave_idx_type id, const Matrix& X, const Material* material, const int32NDArray& nodes)
                :Element3D(id, X, material, nodes) {
                FEM_ASSERT(nodes.numel() == 8);
        }

        virtual const IntegrationRule& GetIntegrationRule(MatrixType eMatType) const {
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

                const double r = rv(0);
                const double s = rv(1);
                const double t = rv(2);

                J(0,0) = (-(X(0,1)*(s+1)*(t+1))/8.0)+(X(0,0)*(s+1)*(t+1))/8.0+(X(0,3)*(1-s)*(t+1))/8.0-(X(0,2)*(1-s)*(t+1))/8.0-(X(0,5)*(s+1)*(1-t))/8.0+(X(0,4)*(s+1)*(1-t))/8.0+(X(0,7)*(1-s)*(1-t))/8.0-(X(0,6)*(1-s)*(1-t))/8.0;
                J(1,0) = (-(X(0,3)*(r+1)*(t+1))/8.0)+(X(0,0)*(r+1)*(t+1))/8.0-(X(0,2)*(1-r)*(t+1))/8.0+(X(0,1)*(1-r)*(t+1))/8.0-(X(0,7)*(r+1)*(1-t))/8.0+(X(0,4)*(r+1)*(1-t))/8.0-(X(0,6)*(1-r)*(1-t))/8.0+(X(0,5)*(1-r)*(1-t))/8.0;
                J(2,0) = (-(X(0,4)*(r+1)*(s+1))/8.0)+(X(0,0)*(r+1)*(s+1))/8.0-(X(0,5)*(1-r)*(s+1))/8.0+(X(0,1)*(1-r)*(s+1))/8.0-(X(0,7)*(r+1)*(1-s))/8.0+(X(0,3)*(r+1)*(1-s))/8.0-(X(0,6)*(1-r)*(1-s))/8.0+(X(0,2)*(1-r)*(1-s))/8.0;
                J(0,1) = (-(X(1,1)*(s+1)*(t+1))/8.0)+(X(1,0)*(s+1)*(t+1))/8.0+(X(1,3)*(1-s)*(t+1))/8.0-(X(1,2)*(1-s)*(t+1))/8.0-(X(1,5)*(s+1)*(1-t))/8.0+(X(1,4)*(s+1)*(1-t))/8.0+(X(1,7)*(1-s)*(1-t))/8.0-(X(1,6)*(1-s)*(1-t))/8.0;
                J(1,1) = (-(X(1,3)*(r+1)*(t+1))/8.0)+(X(1,0)*(r+1)*(t+1))/8.0-(X(1,2)*(1-r)*(t+1))/8.0+(X(1,1)*(1-r)*(t+1))/8.0-(X(1,7)*(r+1)*(1-t))/8.0+(X(1,4)*(r+1)*(1-t))/8.0-(X(1,6)*(1-r)*(1-t))/8.0+(X(1,5)*(1-r)*(1-t))/8.0;
                J(2,1) = (-(X(1,4)*(r+1)*(s+1))/8.0)+(X(1,0)*(r+1)*(s+1))/8.0-(X(1,5)*(1-r)*(s+1))/8.0+(X(1,1)*(1-r)*(s+1))/8.0-(X(1,7)*(r+1)*(1-s))/8.0+(X(1,3)*(r+1)*(1-s))/8.0-(X(1,6)*(1-r)*(1-s))/8.0+(X(1,2)*(1-r)*(1-s))/8.0;
                J(0,2) = (-(X(2,1)*(s+1)*(t+1))/8.0)+(X(2,0)*(s+1)*(t+1))/8.0+(X(2,3)*(1-s)*(t+1))/8.0-(X(2,2)*(1-s)*(t+1))/8.0-(X(2,5)*(s+1)*(1-t))/8.0+(X(2,4)*(s+1)*(1-t))/8.0+(X(2,7)*(1-s)*(1-t))/8.0-(X(2,6)*(1-s)*(1-t))/8.0;
                J(1,2) = (-(X(2,3)*(r+1)*(t+1))/8.0)+(X(2,0)*(r+1)*(t+1))/8.0-(X(2,2)*(1-r)*(t+1))/8.0+(X(2,1)*(1-r)*(t+1))/8.0-(X(2,7)*(r+1)*(1-t))/8.0+(X(2,4)*(r+1)*(1-t))/8.0-(X(2,6)*(1-r)*(1-t))/8.0+(X(2,5)*(1-r)*(1-t))/8.0;
                J(2,2) = (-(X(2,4)*(r+1)*(s+1))/8.0)+(X(2,0)*(r+1)*(s+1))/8.0-(X(2,5)*(1-r)*(s+1))/8.0+(X(2,1)*(1-r)*(s+1))/8.0-(X(2,7)*(r+1)*(1-s))/8.0+(X(2,3)*(r+1)*(1-s))/8.0-(X(2,6)*(1-r)*(1-s))/8.0+(X(2,2)*(1-r)*(1-s))/8.0;

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

                const double r = rv(0);
                const double s = rv(1);
                const double t = rv(2);

                H(0,0) = ((r+1)*(s+1)*(t+1))/8.0;
                H(1,0) = 0;
                H(2,0) = 0;
                H(0,1) = 0;
                H(1,1) = ((r+1)*(s+1)*(t+1))/8.0;
                H(2,1) = 0;
                H(0,2) = 0;
                H(1,2) = 0;
                H(2,2) = ((r+1)*(s+1)*(t+1))/8.0;
                H(0,3) = ((1-r)*(s+1)*(t+1))/8.0;
                H(1,3) = 0;
                H(2,3) = 0;
                H(0,4) = 0;
                H(1,4) = ((1-r)*(s+1)*(t+1))/8.0;
                H(2,4) = 0;
                H(0,5) = 0;
                H(1,5) = 0;
                H(2,5) = ((1-r)*(s+1)*(t+1))/8.0;
                H(0,6) = ((1-r)*(1-s)*(t+1))/8.0;
                H(1,6) = 0;
                H(2,6) = 0;
                H(0,7) = 0;
                H(1,7) = ((1-r)*(1-s)*(t+1))/8.0;
                H(2,7) = 0;
                H(0,8) = 0;
                H(1,8) = 0;
                H(2,8) = ((1-r)*(1-s)*(t+1))/8.0;
                H(0,9) = ((r+1)*(1-s)*(t+1))/8.0;
                H(1,9) = 0;
                H(2,9) = 0;
                H(0,10) = 0;
                H(1,10) = ((r+1)*(1-s)*(t+1))/8.0;
                H(2,10) = 0;
                H(0,11) = 0;
                H(1,11) = 0;
                H(2,11) = ((r+1)*(1-s)*(t+1))/8.0;
                H(0,12) = ((r+1)*(s+1)*(1-t))/8.0;
                H(1,12) = 0;
                H(2,12) = 0;
                H(0,13) = 0;
                H(1,13) = ((r+1)*(s+1)*(1-t))/8.0;
                H(2,13) = 0;
                H(0,14) = 0;
                H(1,14) = 0;
                H(2,14) = ((r+1)*(s+1)*(1-t))/8.0;
                H(0,15) = ((1-r)*(s+1)*(1-t))/8.0;
                H(1,15) = 0;
                H(2,15) = 0;
                H(0,16) = 0;
                H(1,16) = ((1-r)*(s+1)*(1-t))/8.0;
                H(2,16) = 0;
                H(0,17) = 0;
                H(1,17) = 0;
                H(2,17) = ((1-r)*(s+1)*(1-t))/8.0;
                H(0,18) = ((1-r)*(1-s)*(1-t))/8.0;
                H(1,18) = 0;
                H(2,18) = 0;
                H(0,19) = 0;
                H(1,19) = ((1-r)*(1-s)*(1-t))/8.0;
                H(2,19) = 0;
                H(0,20) = 0;
                H(1,20) = 0;
                H(2,20) = ((1-r)*(1-s)*(1-t))/8.0;
                H(0,21) = ((r+1)*(1-s)*(1-t))/8.0;
                H(1,21) = 0;
                H(2,21) = 0;
                H(0,22) = 0;
                H(1,22) = ((r+1)*(1-s)*(1-t))/8.0;
                H(2,22) = 0;
                H(0,23) = 0;
                H(1,23) = 0;
                H(2,23) = ((r+1)*(1-s)*(1-t))/8.0;
        }

        virtual void StrainMatrix(const ColumnVector& rv, const Matrix& invJ, Matrix& B) const {
                FEM_ASSERT(rv.numel() == 3);
                FEM_ASSERT(invJ.rows() == 3);
                FEM_ASSERT(invJ.columns() == 3);
                FEM_ASSERT(B.rows() == 6);
                FEM_ASSERT(B.columns() == 24);

                const double r = rv(0);
                const double s = rv(1);
                const double t = rv(2);

                B(0,0) = (invJ(0,0)*(s+1)*(t+1))/8.0+(invJ(0,1)*(r+1)*(t+1))/8.0+(invJ(0,2)*(r+1)*(s+1))/8.0;
                B(1,0) = 0;
                B(2,0) = 0;
                B(3,0) = (invJ(1,0)*(s+1)*(t+1))/8.0+(invJ(1,1)*(r+1)*(t+1))/8.0+(invJ(1,2)*(r+1)*(s+1))/8.0;
                B(4,0) = 0;
                B(5,0) = (invJ(2,0)*(s+1)*(t+1))/8.0+(invJ(2,1)*(r+1)*(t+1))/8.0+(invJ(2,2)*(r+1)*(s+1))/8.0;
                B(0,1) = 0;
                B(1,1) = (invJ(1,0)*(s+1)*(t+1))/8.0+(invJ(1,1)*(r+1)*(t+1))/8.0+(invJ(1,2)*(r+1)*(s+1))/8.0;
                B(2,1) = 0;
                B(3,1) = (invJ(0,0)*(s+1)*(t+1))/8.0+(invJ(0,1)*(r+1)*(t+1))/8.0+(invJ(0,2)*(r+1)*(s+1))/8.0;
                B(4,1) = (invJ(2,0)*(s+1)*(t+1))/8.0+(invJ(2,1)*(r+1)*(t+1))/8.0+(invJ(2,2)*(r+1)*(s+1))/8.0;
                B(5,1) = 0;
                B(0,2) = 0;
                B(1,2) = 0;
                B(2,2) = (invJ(2,0)*(s+1)*(t+1))/8.0+(invJ(2,1)*(r+1)*(t+1))/8.0+(invJ(2,2)*(r+1)*(s+1))/8.0;
                B(3,2) = 0;
                B(4,2) = (invJ(1,0)*(s+1)*(t+1))/8.0+(invJ(1,1)*(r+1)*(t+1))/8.0+(invJ(1,2)*(r+1)*(s+1))/8.0;
                B(5,2) = (invJ(0,0)*(s+1)*(t+1))/8.0+(invJ(0,1)*(r+1)*(t+1))/8.0+(invJ(0,2)*(r+1)*(s+1))/8.0;
                B(0,3) = (-(invJ(0,0)*(s+1)*(t+1))/8.0)+(invJ(0,1)*(1-r)*(t+1))/8.0+(invJ(0,2)*(1-r)*(s+1))/8.0;
                B(1,3) = 0;
                B(2,3) = 0;
                B(3,3) = (-(invJ(1,0)*(s+1)*(t+1))/8.0)+(invJ(1,1)*(1-r)*(t+1))/8.0+(invJ(1,2)*(1-r)*(s+1))/8.0;
                B(4,3) = 0;
                B(5,3) = (-(invJ(2,0)*(s+1)*(t+1))/8.0)+(invJ(2,1)*(1-r)*(t+1))/8.0+(invJ(2,2)*(1-r)*(s+1))/8.0;
                B(0,4) = 0;
                B(1,4) = (-(invJ(1,0)*(s+1)*(t+1))/8.0)+(invJ(1,1)*(1-r)*(t+1))/8.0+(invJ(1,2)*(1-r)*(s+1))/8.0;
                B(2,4) = 0;
                B(3,4) = (-(invJ(0,0)*(s+1)*(t+1))/8.0)+(invJ(0,1)*(1-r)*(t+1))/8.0+(invJ(0,2)*(1-r)*(s+1))/8.0;
                B(4,4) = (-(invJ(2,0)*(s+1)*(t+1))/8.0)+(invJ(2,1)*(1-r)*(t+1))/8.0+(invJ(2,2)*(1-r)*(s+1))/8.0;
                B(5,4) = 0;
                B(0,5) = 0;
                B(1,5) = 0;
                B(2,5) = (-(invJ(2,0)*(s+1)*(t+1))/8.0)+(invJ(2,1)*(1-r)*(t+1))/8.0+(invJ(2,2)*(1-r)*(s+1))/8.0;
                B(3,5) = 0;
                B(4,5) = (-(invJ(1,0)*(s+1)*(t+1))/8.0)+(invJ(1,1)*(1-r)*(t+1))/8.0+(invJ(1,2)*(1-r)*(s+1))/8.0;
                B(5,5) = (-(invJ(0,0)*(s+1)*(t+1))/8.0)+(invJ(0,1)*(1-r)*(t+1))/8.0+(invJ(0,2)*(1-r)*(s+1))/8.0;
                B(0,6) = (-(invJ(0,0)*(1-s)*(t+1))/8.0)-(invJ(0,1)*(1-r)*(t+1))/8.0+(invJ(0,2)*(1-r)*(1-s))/8.0;
                B(1,6) = 0;
                B(2,6) = 0;
                B(3,6) = (-(invJ(1,0)*(1-s)*(t+1))/8.0)-(invJ(1,1)*(1-r)*(t+1))/8.0+(invJ(1,2)*(1-r)*(1-s))/8.0;
                B(4,6) = 0;
                B(5,6) = (-(invJ(2,0)*(1-s)*(t+1))/8.0)-(invJ(2,1)*(1-r)*(t+1))/8.0+(invJ(2,2)*(1-r)*(1-s))/8.0;
                B(0,7) = 0;
                B(1,7) = (-(invJ(1,0)*(1-s)*(t+1))/8.0)-(invJ(1,1)*(1-r)*(t+1))/8.0+(invJ(1,2)*(1-r)*(1-s))/8.0;
                B(2,7) = 0;
                B(3,7) = (-(invJ(0,0)*(1-s)*(t+1))/8.0)-(invJ(0,1)*(1-r)*(t+1))/8.0+(invJ(0,2)*(1-r)*(1-s))/8.0;
                B(4,7) = (-(invJ(2,0)*(1-s)*(t+1))/8.0)-(invJ(2,1)*(1-r)*(t+1))/8.0+(invJ(2,2)*(1-r)*(1-s))/8.0;
                B(5,7) = 0;
                B(0,8) = 0;
                B(1,8) = 0;
                B(2,8) = (-(invJ(2,0)*(1-s)*(t+1))/8.0)-(invJ(2,1)*(1-r)*(t+1))/8.0+(invJ(2,2)*(1-r)*(1-s))/8.0;
                B(3,8) = 0;
                B(4,8) = (-(invJ(1,0)*(1-s)*(t+1))/8.0)-(invJ(1,1)*(1-r)*(t+1))/8.0+(invJ(1,2)*(1-r)*(1-s))/8.0;
                B(5,8) = (-(invJ(0,0)*(1-s)*(t+1))/8.0)-(invJ(0,1)*(1-r)*(t+1))/8.0+(invJ(0,2)*(1-r)*(1-s))/8.0;
                B(0,9) = (invJ(0,0)*(1-s)*(t+1))/8.0-(invJ(0,1)*(r+1)*(t+1))/8.0+(invJ(0,2)*(r+1)*(1-s))/8.0;
                B(1,9) = 0;
                B(2,9) = 0;
                B(3,9) = (invJ(1,0)*(1-s)*(t+1))/8.0-(invJ(1,1)*(r+1)*(t+1))/8.0+(invJ(1,2)*(r+1)*(1-s))/8.0;
                B(4,9) = 0;
                B(5,9) = (invJ(2,0)*(1-s)*(t+1))/8.0-(invJ(2,1)*(r+1)*(t+1))/8.0+(invJ(2,2)*(r+1)*(1-s))/8.0;
                B(0,10) = 0;
                B(1,10) = (invJ(1,0)*(1-s)*(t+1))/8.0-(invJ(1,1)*(r+1)*(t+1))/8.0+(invJ(1,2)*(r+1)*(1-s))/8.0;
                B(2,10) = 0;
                B(3,10) = (invJ(0,0)*(1-s)*(t+1))/8.0-(invJ(0,1)*(r+1)*(t+1))/8.0+(invJ(0,2)*(r+1)*(1-s))/8.0;
                B(4,10) = (invJ(2,0)*(1-s)*(t+1))/8.0-(invJ(2,1)*(r+1)*(t+1))/8.0+(invJ(2,2)*(r+1)*(1-s))/8.0;
                B(5,10) = 0;
                B(0,11) = 0;
                B(1,11) = 0;
                B(2,11) = (invJ(2,0)*(1-s)*(t+1))/8.0-(invJ(2,1)*(r+1)*(t+1))/8.0+(invJ(2,2)*(r+1)*(1-s))/8.0;
                B(3,11) = 0;
                B(4,11) = (invJ(1,0)*(1-s)*(t+1))/8.0-(invJ(1,1)*(r+1)*(t+1))/8.0+(invJ(1,2)*(r+1)*(1-s))/8.0;
                B(5,11) = (invJ(0,0)*(1-s)*(t+1))/8.0-(invJ(0,1)*(r+1)*(t+1))/8.0+(invJ(0,2)*(r+1)*(1-s))/8.0;
                B(0,12) = (invJ(0,0)*(s+1)*(1-t))/8.0+(invJ(0,1)*(r+1)*(1-t))/8.0-(invJ(0,2)*(r+1)*(s+1))/8.0;
                B(1,12) = 0;
                B(2,12) = 0;
                B(3,12) = (invJ(1,0)*(s+1)*(1-t))/8.0+(invJ(1,1)*(r+1)*(1-t))/8.0-(invJ(1,2)*(r+1)*(s+1))/8.0;
                B(4,12) = 0;
                B(5,12) = (invJ(2,0)*(s+1)*(1-t))/8.0+(invJ(2,1)*(r+1)*(1-t))/8.0-(invJ(2,2)*(r+1)*(s+1))/8.0;
                B(0,13) = 0;
                B(1,13) = (invJ(1,0)*(s+1)*(1-t))/8.0+(invJ(1,1)*(r+1)*(1-t))/8.0-(invJ(1,2)*(r+1)*(s+1))/8.0;
                B(2,13) = 0;
                B(3,13) = (invJ(0,0)*(s+1)*(1-t))/8.0+(invJ(0,1)*(r+1)*(1-t))/8.0-(invJ(0,2)*(r+1)*(s+1))/8.0;
                B(4,13) = (invJ(2,0)*(s+1)*(1-t))/8.0+(invJ(2,1)*(r+1)*(1-t))/8.0-(invJ(2,2)*(r+1)*(s+1))/8.0;
                B(5,13) = 0;
                B(0,14) = 0;
                B(1,14) = 0;
                B(2,14) = (invJ(2,0)*(s+1)*(1-t))/8.0+(invJ(2,1)*(r+1)*(1-t))/8.0-(invJ(2,2)*(r+1)*(s+1))/8.0;
                B(3,14) = 0;
                B(4,14) = (invJ(1,0)*(s+1)*(1-t))/8.0+(invJ(1,1)*(r+1)*(1-t))/8.0-(invJ(1,2)*(r+1)*(s+1))/8.0;
                B(5,14) = (invJ(0,0)*(s+1)*(1-t))/8.0+(invJ(0,1)*(r+1)*(1-t))/8.0-(invJ(0,2)*(r+1)*(s+1))/8.0;
                B(0,15) = (-(invJ(0,0)*(s+1)*(1-t))/8.0)+(invJ(0,1)*(1-r)*(1-t))/8.0-(invJ(0,2)*(1-r)*(s+1))/8.0;
                B(1,15) = 0;
                B(2,15) = 0;
                B(3,15) = (-(invJ(1,0)*(s+1)*(1-t))/8.0)+(invJ(1,1)*(1-r)*(1-t))/8.0-(invJ(1,2)*(1-r)*(s+1))/8.0;
                B(4,15) = 0;
                B(5,15) = (-(invJ(2,0)*(s+1)*(1-t))/8.0)+(invJ(2,1)*(1-r)*(1-t))/8.0-(invJ(2,2)*(1-r)*(s+1))/8.0;
                B(0,16) = 0;
                B(1,16) = (-(invJ(1,0)*(s+1)*(1-t))/8.0)+(invJ(1,1)*(1-r)*(1-t))/8.0-(invJ(1,2)*(1-r)*(s+1))/8.0;
                B(2,16) = 0;
                B(3,16) = (-(invJ(0,0)*(s+1)*(1-t))/8.0)+(invJ(0,1)*(1-r)*(1-t))/8.0-(invJ(0,2)*(1-r)*(s+1))/8.0;
                B(4,16) = (-(invJ(2,0)*(s+1)*(1-t))/8.0)+(invJ(2,1)*(1-r)*(1-t))/8.0-(invJ(2,2)*(1-r)*(s+1))/8.0;
                B(5,16) = 0;
                B(0,17) = 0;
                B(1,17) = 0;
                B(2,17) = (-(invJ(2,0)*(s+1)*(1-t))/8.0)+(invJ(2,1)*(1-r)*(1-t))/8.0-(invJ(2,2)*(1-r)*(s+1))/8.0;
                B(3,17) = 0;
                B(4,17) = (-(invJ(1,0)*(s+1)*(1-t))/8.0)+(invJ(1,1)*(1-r)*(1-t))/8.0-(invJ(1,2)*(1-r)*(s+1))/8.0;
                B(5,17) = (-(invJ(0,0)*(s+1)*(1-t))/8.0)+(invJ(0,1)*(1-r)*(1-t))/8.0-(invJ(0,2)*(1-r)*(s+1))/8.0;
                B(0,18) = (-(invJ(0,0)*(1-s)*(1-t))/8.0)-(invJ(0,1)*(1-r)*(1-t))/8.0-(invJ(0,2)*(1-r)*(1-s))/8.0;
                B(1,18) = 0;
                B(2,18) = 0;
                B(3,18) = (-(invJ(1,0)*(1-s)*(1-t))/8.0)-(invJ(1,1)*(1-r)*(1-t))/8.0-(invJ(1,2)*(1-r)*(1-s))/8.0;
                B(4,18) = 0;
                B(5,18) = (-(invJ(2,0)*(1-s)*(1-t))/8.0)-(invJ(2,1)*(1-r)*(1-t))/8.0-(invJ(2,2)*(1-r)*(1-s))/8.0;
                B(0,19) = 0;
                B(1,19) = (-(invJ(1,0)*(1-s)*(1-t))/8.0)-(invJ(1,1)*(1-r)*(1-t))/8.0-(invJ(1,2)*(1-r)*(1-s))/8.0;
                B(2,19) = 0;
                B(3,19) = (-(invJ(0,0)*(1-s)*(1-t))/8.0)-(invJ(0,1)*(1-r)*(1-t))/8.0-(invJ(0,2)*(1-r)*(1-s))/8.0;
                B(4,19) = (-(invJ(2,0)*(1-s)*(1-t))/8.0)-(invJ(2,1)*(1-r)*(1-t))/8.0-(invJ(2,2)*(1-r)*(1-s))/8.0;
                B(5,19) = 0;
                B(0,20) = 0;
                B(1,20) = 0;
                B(2,20) = (-(invJ(2,0)*(1-s)*(1-t))/8.0)-(invJ(2,1)*(1-r)*(1-t))/8.0-(invJ(2,2)*(1-r)*(1-s))/8.0;
                B(3,20) = 0;
                B(4,20) = (-(invJ(1,0)*(1-s)*(1-t))/8.0)-(invJ(1,1)*(1-r)*(1-t))/8.0-(invJ(1,2)*(1-r)*(1-s))/8.0;
                B(5,20) = (-(invJ(0,0)*(1-s)*(1-t))/8.0)-(invJ(0,1)*(1-r)*(1-t))/8.0-(invJ(0,2)*(1-r)*(1-s))/8.0;
                B(0,21) = (invJ(0,0)*(1-s)*(1-t))/8.0-(invJ(0,1)*(r+1)*(1-t))/8.0-(invJ(0,2)*(r+1)*(1-s))/8.0;
                B(1,21) = 0;
                B(2,21) = 0;
                B(3,21) = (invJ(1,0)*(1-s)*(1-t))/8.0-(invJ(1,1)*(r+1)*(1-t))/8.0-(invJ(1,2)*(r+1)*(1-s))/8.0;
                B(4,21) = 0;
                B(5,21) = (invJ(2,0)*(1-s)*(1-t))/8.0-(invJ(2,1)*(r+1)*(1-t))/8.0-(invJ(2,2)*(r+1)*(1-s))/8.0;
                B(0,22) = 0;
                B(1,22) = (invJ(1,0)*(1-s)*(1-t))/8.0-(invJ(1,1)*(r+1)*(1-t))/8.0-(invJ(1,2)*(r+1)*(1-s))/8.0;
                B(2,22) = 0;
                B(3,22) = (invJ(0,0)*(1-s)*(1-t))/8.0-(invJ(0,1)*(r+1)*(1-t))/8.0-(invJ(0,2)*(r+1)*(1-s))/8.0;
                B(4,22) = (invJ(2,0)*(1-s)*(1-t))/8.0-(invJ(2,1)*(r+1)*(1-t))/8.0-(invJ(2,2)*(r+1)*(1-s))/8.0;
                B(5,22) = 0;
                B(0,23) = 0;
                B(1,23) = 0;
                B(2,23) = (invJ(2,0)*(1-s)*(1-t))/8.0-(invJ(2,1)*(r+1)*(1-t))/8.0-(invJ(2,2)*(r+1)*(1-s))/8.0;
                B(3,23) = 0;
                B(4,23) = (invJ(1,0)*(1-s)*(1-t))/8.0-(invJ(1,1)*(r+1)*(1-t))/8.0-(invJ(1,2)*(r+1)*(1-s))/8.0;
                B(5,23) = (invJ(0,0)*(1-s)*(1-t))/8.0-(invJ(0,1)*(r+1)*(1-t))/8.0-(invJ(0,2)*(r+1)*(1-s))/8.0;

#ifdef DEBUG
                for (octave_idx_type i = 0; i < B.rows(); ++i) {
                        for (octave_idx_type j = 0; j < B.columns(); ++j) {
                                FEM_ASSERT(std::isfinite(B(i, j)));
                        }
                }
#endif
        }

        virtual Matrix InterpGaussToNodal(MatrixType eMatType, const Matrix& taug) const {
                const IntegrationRule& oIntegRule = GetIntegrationRule(eMatType);
                const octave_idx_type iNumGauss = oIntegRule.iGetNumEvalPoints();
                const octave_idx_type iNumDir = oIntegRule.iGetNumDirections();
                const octave_idx_type iNumNodes = nodes.numel();
                ColumnVector rv(iNumDir);

                Matrix H(iNumGauss, iNumNodes);

                for (octave_idx_type i = 0; i < iNumGauss; ++i) {
                        for (octave_idx_type j = 0; j < iNumDir; ++j) {
                                rv(j) = oIntegRule.dGetPosition(i, j);
                        }

                        StressInterpMatrix(rv, H, i);
                }

                return H.solve(taug);
        }

private:
        void StressInterpMatrix(const ColumnVector& rv, Matrix& Hs, octave_idx_type irow) const {
                FEM_ASSERT(rv.numel() == 3);
                FEM_ASSERT(Hs.columns() == 8);
                FEM_ASSERT(irow >= 0);
                FEM_ASSERT(irow < Hs.rows());

                const double r = rv(0);
                const double s = rv(1);
                const double t = rv(2);

                Hs(irow, 0) = ((r+1)*(s+1)*(t+1))/8.0;
                Hs(irow, 1) = ((1-r)*(s+1)*(t+1))/8.0;
                Hs(irow, 2) = ((1-r)*(1-s)*(t+1))/8.0;
                Hs(irow, 3) = ((r+1)*(1-s)*(t+1))/8.0;
                Hs(irow, 4) = ((r+1)*(s+1)*(1-t))/8.0;
                Hs(irow, 5) = ((1-r)*(s+1)*(1-t))/8.0;
                Hs(irow, 6) = ((1-r)*(1-s)*(1-t))/8.0;
                Hs(irow, 7) = ((r+1)*(1-s)*(1-t))/8.0;
        }
};

class Tet10: public Element3D
{
        static constexpr double gamma = 1. / 6.;

public:
        Tet10(octave_idx_type id, const Matrix& X, const Material* material, const int32NDArray& nodes)
                :Element3D(id, X, material, nodes) {
                FEM_ASSERT(nodes.numel() == 10);
        }

        virtual const IntegrationRule& GetIntegrationRule(MatrixType eMatType) const {
                static IntegrationRule oIntegStiff, oIntegMass, oIntegMassDiag;

                switch (eMatType) {
                case MAT_STIFFNESS:
                case MAT_STIFFNESS_SYM:
                case MAT_STIFFNESS_SYM_L:
                case VEC_STRESS_CAUCH:
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

                const double Zeta1 = rv(0);
                const double Zeta2 = rv(1);
                const double Zeta3 = rv(2);
                const double Zeta4 = rv(3);

                FEM_ASSERT(Zeta1 + Zeta2 + Zeta3 + Zeta4 == 1);

                J(0,0) = 1;
                J(1,0) = 4*X(0,7)*Zeta4+4*X(0,6)*Zeta3+4*X(0,4)*Zeta2+X(0,0)*(4*Zeta1-1);
                J(2,0) = 4*X(1,7)*Zeta4+4*X(1,6)*Zeta3+4*X(1,4)*Zeta2+X(1,0)*(4*Zeta1-1);
                J(3,0) = 4*X(2,7)*Zeta4+4*X(2,6)*Zeta3+4*X(2,4)*Zeta2+X(2,0)*(4*Zeta1-1);
                J(0,1) = 1;
                J(1,1) = 4*X(0,8)*Zeta4+4*X(0,5)*Zeta3+X(0,1)*(4*Zeta2-1)+4*X(0,4)*Zeta1;
                J(2,1) = 4*X(1,8)*Zeta4+4*X(1,5)*Zeta3+X(1,1)*(4*Zeta2-1)+4*X(1,4)*Zeta1;
                J(3,1) = 4*X(2,8)*Zeta4+4*X(2,5)*Zeta3+X(2,1)*(4*Zeta2-1)+4*X(2,4)*Zeta1;
                J(0,2) = 1;
                J(1,2) = 4*X(0,9)*Zeta4+X(0,2)*(4*Zeta3-1)+4*X(0,5)*Zeta2+4*X(0,6)*Zeta1;
                J(2,2) = 4*X(1,9)*Zeta4+X(1,2)*(4*Zeta3-1)+4*X(1,5)*Zeta2+4*X(1,6)*Zeta1;
                J(3,2) = 4*X(2,9)*Zeta4+X(2,2)*(4*Zeta3-1)+4*X(2,5)*Zeta2+4*X(2,6)*Zeta1;
                J(0,3) = 1;
                J(1,3) = X(0,3)*(4*Zeta4-1)+4*X(0,9)*Zeta3+4*X(0,8)*Zeta2+4*X(0,7)*Zeta1;
                J(2,3) = X(1,3)*(4*Zeta4-1)+4*X(1,9)*Zeta3+4*X(1,8)*Zeta2+4*X(1,7)*Zeta1;
                J(3,3) = X(2,3)*(4*Zeta4-1)+4*X(2,9)*Zeta3+4*X(2,8)*Zeta2+4*X(2,7)*Zeta1;

                return J.determinant() * gamma;
        }

        virtual void DispInterpMatrix(const ColumnVector& rv, Matrix& H) const {
                FEM_ASSERT(H.rows() == 3);
                FEM_ASSERT(H.columns() == 30);
                FEM_ASSERT(rv.numel() == 4);

                const double Zeta1 = rv(0);
                const double Zeta2 = rv(1);
                const double Zeta3 = rv(2);
                const double Zeta4 = rv(3);

                H(0,0) = Zeta1*(2*Zeta1-1);
                H(1,0) = 0;
                H(2,0) = 0;
                H(0,1) = 0;
                H(1,1) = Zeta1*(2*Zeta1-1);
                H(2,1) = 0;
                H(0,2) = 0;
                H(1,2) = 0;
                H(2,2) = Zeta1*(2*Zeta1-1);
                H(0,3) = Zeta2*(2*Zeta2-1);
                H(1,3) = 0;
                H(2,3) = 0;
                H(0,4) = 0;
                H(1,4) = Zeta2*(2*Zeta2-1);
                H(2,4) = 0;
                H(0,5) = 0;
                H(1,5) = 0;
                H(2,5) = Zeta2*(2*Zeta2-1);
                H(0,6) = Zeta3*(2*Zeta3-1);
                H(1,6) = 0;
                H(2,6) = 0;
                H(0,7) = 0;
                H(1,7) = Zeta3*(2*Zeta3-1);
                H(2,7) = 0;
                H(0,8) = 0;
                H(1,8) = 0;
                H(2,8) = Zeta3*(2*Zeta3-1);
                H(0,9) = Zeta4*(2*Zeta4-1);
                H(1,9) = 0;
                H(2,9) = 0;
                H(0,10) = 0;
                H(1,10) = Zeta4*(2*Zeta4-1);
                H(2,10) = 0;
                H(0,11) = 0;
                H(1,11) = 0;
                H(2,11) = Zeta4*(2*Zeta4-1);
                H(0,12) = 4*Zeta1*Zeta2;
                H(1,12) = 0;
                H(2,12) = 0;
                H(0,13) = 0;
                H(1,13) = 4*Zeta1*Zeta2;
                H(2,13) = 0;
                H(0,14) = 0;
                H(1,14) = 0;
                H(2,14) = 4*Zeta1*Zeta2;
                H(0,15) = 4*Zeta2*Zeta3;
                H(1,15) = 0;
                H(2,15) = 0;
                H(0,16) = 0;
                H(1,16) = 4*Zeta2*Zeta3;
                H(2,16) = 0;
                H(0,17) = 0;
                H(1,17) = 0;
                H(2,17) = 4*Zeta2*Zeta3;
                H(0,18) = 4*Zeta1*Zeta3;
                H(1,18) = 0;
                H(2,18) = 0;
                H(0,19) = 0;
                H(1,19) = 4*Zeta1*Zeta3;
                H(2,19) = 0;
                H(0,20) = 0;
                H(1,20) = 0;
                H(2,20) = 4*Zeta1*Zeta3;
                H(0,21) = 4*Zeta1*Zeta4;
                H(1,21) = 0;
                H(2,21) = 0;
                H(0,22) = 0;
                H(1,22) = 4*Zeta1*Zeta4;
                H(2,22) = 0;
                H(0,23) = 0;
                H(1,23) = 0;
                H(2,23) = 4*Zeta1*Zeta4;
                H(0,24) = 4*Zeta2*Zeta4;
                H(1,24) = 0;
                H(2,24) = 0;
                H(0,25) = 0;
                H(1,25) = 4*Zeta2*Zeta4;
                H(2,25) = 0;
                H(0,26) = 0;
                H(1,26) = 0;
                H(2,26) = 4*Zeta2*Zeta4;
                H(0,27) = 4*Zeta3*Zeta4;
                H(1,27) = 0;
                H(2,27) = 0;
                H(0,28) = 0;
                H(1,28) = 4*Zeta3*Zeta4;
                H(2,28) = 0;
                H(0,29) = 0;
                H(1,29) = 0;
                H(2,29) = 4*Zeta3*Zeta4;
        }

        virtual void StrainMatrix(const ColumnVector& rv, const Matrix& invJ, Matrix& B) const {
                FEM_ASSERT(invJ.rows() == 4);
                FEM_ASSERT(invJ.columns() == 4);
                FEM_ASSERT(B.rows() == 6);
                FEM_ASSERT(B.columns() == 30);
                FEM_ASSERT(rv.numel() == 4);

                const double Zeta1 = rv(0);
                const double Zeta2 = rv(1);
                const double Zeta3 = rv(2);
                const double Zeta4 = rv(3);

                B(0,0) = invJ(0,1)*(4*Zeta1-1);
                B(1,0) = 0;
                B(2,0) = 0;
                B(3,0) = invJ(0,2)*(4*Zeta1-1);
                B(4,0) = 0;
                B(5,0) = invJ(0,3)*(4*Zeta1-1);
                B(0,1) = 0;
                B(1,1) = invJ(0,2)*(4*Zeta1-1);
                B(2,1) = 0;
                B(3,1) = invJ(0,1)*(4*Zeta1-1);
                B(4,1) = invJ(0,3)*(4*Zeta1-1);
                B(5,1) = 0;
                B(0,2) = 0;
                B(1,2) = 0;
                B(2,2) = invJ(0,3)*(4*Zeta1-1);
                B(3,2) = 0;
                B(4,2) = invJ(0,2)*(4*Zeta1-1);
                B(5,2) = invJ(0,1)*(4*Zeta1-1);
                B(0,3) = invJ(1,1)*(4*Zeta2-1);
                B(1,3) = 0;
                B(2,3) = 0;
                B(3,3) = invJ(1,2)*(4*Zeta2-1);
                B(4,3) = 0;
                B(5,3) = invJ(1,3)*(4*Zeta2-1);
                B(0,4) = 0;
                B(1,4) = invJ(1,2)*(4*Zeta2-1);
                B(2,4) = 0;
                B(3,4) = invJ(1,1)*(4*Zeta2-1);
                B(4,4) = invJ(1,3)*(4*Zeta2-1);
                B(5,4) = 0;
                B(0,5) = 0;
                B(1,5) = 0;
                B(2,5) = invJ(1,3)*(4*Zeta2-1);
                B(3,5) = 0;
                B(4,5) = invJ(1,2)*(4*Zeta2-1);
                B(5,5) = invJ(1,1)*(4*Zeta2-1);
                B(0,6) = invJ(2,1)*(4*Zeta3-1);
                B(1,6) = 0;
                B(2,6) = 0;
                B(3,6) = invJ(2,2)*(4*Zeta3-1);
                B(4,6) = 0;
                B(5,6) = invJ(2,3)*(4*Zeta3-1);
                B(0,7) = 0;
                B(1,7) = invJ(2,2)*(4*Zeta3-1);
                B(2,7) = 0;
                B(3,7) = invJ(2,1)*(4*Zeta3-1);
                B(4,7) = invJ(2,3)*(4*Zeta3-1);
                B(5,7) = 0;
                B(0,8) = 0;
                B(1,8) = 0;
                B(2,8) = invJ(2,3)*(4*Zeta3-1);
                B(3,8) = 0;
                B(4,8) = invJ(2,2)*(4*Zeta3-1);
                B(5,8) = invJ(2,1)*(4*Zeta3-1);
                B(0,9) = invJ(3,1)*(4*Zeta4-1);
                B(1,9) = 0;
                B(2,9) = 0;
                B(3,9) = invJ(3,2)*(4*Zeta4-1);
                B(4,9) = 0;
                B(5,9) = invJ(3,3)*(4*Zeta4-1);
                B(0,10) = 0;
                B(1,10) = invJ(3,2)*(4*Zeta4-1);
                B(2,10) = 0;
                B(3,10) = invJ(3,1)*(4*Zeta4-1);
                B(4,10) = invJ(3,3)*(4*Zeta4-1);
                B(5,10) = 0;
                B(0,11) = 0;
                B(1,11) = 0;
                B(2,11) = invJ(3,3)*(4*Zeta4-1);
                B(3,11) = 0;
                B(4,11) = invJ(3,2)*(4*Zeta4-1);
                B(5,11) = invJ(3,1)*(4*Zeta4-1);
                B(0,12) = 4*invJ(0,1)*Zeta2+4*invJ(1,1)*Zeta1;
                B(1,12) = 0;
                B(2,12) = 0;
                B(3,12) = 4*invJ(0,2)*Zeta2+4*invJ(1,2)*Zeta1;
                B(4,12) = 0;
                B(5,12) = 4*invJ(0,3)*Zeta2+4*invJ(1,3)*Zeta1;
                B(0,13) = 0;
                B(1,13) = 4*invJ(0,2)*Zeta2+4*invJ(1,2)*Zeta1;
                B(2,13) = 0;
                B(3,13) = 4*invJ(0,1)*Zeta2+4*invJ(1,1)*Zeta1;
                B(4,13) = 4*invJ(0,3)*Zeta2+4*invJ(1,3)*Zeta1;
                B(5,13) = 0;
                B(0,14) = 0;
                B(1,14) = 0;
                B(2,14) = 4*invJ(0,3)*Zeta2+4*invJ(1,3)*Zeta1;
                B(3,14) = 0;
                B(4,14) = 4*invJ(0,2)*Zeta2+4*invJ(1,2)*Zeta1;
                B(5,14) = 4*invJ(0,1)*Zeta2+4*invJ(1,1)*Zeta1;
                B(0,15) = 4*invJ(1,1)*Zeta3+4*invJ(2,1)*Zeta2;
                B(1,15) = 0;
                B(2,15) = 0;
                B(3,15) = 4*invJ(1,2)*Zeta3+4*invJ(2,2)*Zeta2;
                B(4,15) = 0;
                B(5,15) = 4*invJ(1,3)*Zeta3+4*invJ(2,3)*Zeta2;
                B(0,16) = 0;
                B(1,16) = 4*invJ(1,2)*Zeta3+4*invJ(2,2)*Zeta2;
                B(2,16) = 0;
                B(3,16) = 4*invJ(1,1)*Zeta3+4*invJ(2,1)*Zeta2;
                B(4,16) = 4*invJ(1,3)*Zeta3+4*invJ(2,3)*Zeta2;
                B(5,16) = 0;
                B(0,17) = 0;
                B(1,17) = 0;
                B(2,17) = 4*invJ(1,3)*Zeta3+4*invJ(2,3)*Zeta2;
                B(3,17) = 0;
                B(4,17) = 4*invJ(1,2)*Zeta3+4*invJ(2,2)*Zeta2;
                B(5,17) = 4*invJ(1,1)*Zeta3+4*invJ(2,1)*Zeta2;
                B(0,18) = 4*invJ(0,1)*Zeta3+4*invJ(2,1)*Zeta1;
                B(1,18) = 0;
                B(2,18) = 0;
                B(3,18) = 4*invJ(0,2)*Zeta3+4*invJ(2,2)*Zeta1;
                B(4,18) = 0;
                B(5,18) = 4*invJ(0,3)*Zeta3+4*invJ(2,3)*Zeta1;
                B(0,19) = 0;
                B(1,19) = 4*invJ(0,2)*Zeta3+4*invJ(2,2)*Zeta1;
                B(2,19) = 0;
                B(3,19) = 4*invJ(0,1)*Zeta3+4*invJ(2,1)*Zeta1;
                B(4,19) = 4*invJ(0,3)*Zeta3+4*invJ(2,3)*Zeta1;
                B(5,19) = 0;
                B(0,20) = 0;
                B(1,20) = 0;
                B(2,20) = 4*invJ(0,3)*Zeta3+4*invJ(2,3)*Zeta1;
                B(3,20) = 0;
                B(4,20) = 4*invJ(0,2)*Zeta3+4*invJ(2,2)*Zeta1;
                B(5,20) = 4*invJ(0,1)*Zeta3+4*invJ(2,1)*Zeta1;
                B(0,21) = 4*invJ(0,1)*Zeta4+4*invJ(3,1)*Zeta1;
                B(1,21) = 0;
                B(2,21) = 0;
                B(3,21) = 4*invJ(0,2)*Zeta4+4*invJ(3,2)*Zeta1;
                B(4,21) = 0;
                B(5,21) = 4*invJ(0,3)*Zeta4+4*invJ(3,3)*Zeta1;
                B(0,22) = 0;
                B(1,22) = 4*invJ(0,2)*Zeta4+4*invJ(3,2)*Zeta1;
                B(2,22) = 0;
                B(3,22) = 4*invJ(0,1)*Zeta4+4*invJ(3,1)*Zeta1;
                B(4,22) = 4*invJ(0,3)*Zeta4+4*invJ(3,3)*Zeta1;
                B(5,22) = 0;
                B(0,23) = 0;
                B(1,23) = 0;
                B(2,23) = 4*invJ(0,3)*Zeta4+4*invJ(3,3)*Zeta1;
                B(3,23) = 0;
                B(4,23) = 4*invJ(0,2)*Zeta4+4*invJ(3,2)*Zeta1;
                B(5,23) = 4*invJ(0,1)*Zeta4+4*invJ(3,1)*Zeta1;
                B(0,24) = 4*invJ(1,1)*Zeta4+4*invJ(3,1)*Zeta2;
                B(1,24) = 0;
                B(2,24) = 0;
                B(3,24) = 4*invJ(1,2)*Zeta4+4*invJ(3,2)*Zeta2;
                B(4,24) = 0;
                B(5,24) = 4*invJ(1,3)*Zeta4+4*invJ(3,3)*Zeta2;
                B(0,25) = 0;
                B(1,25) = 4*invJ(1,2)*Zeta4+4*invJ(3,2)*Zeta2;
                B(2,25) = 0;
                B(3,25) = 4*invJ(1,1)*Zeta4+4*invJ(3,1)*Zeta2;
                B(4,25) = 4*invJ(1,3)*Zeta4+4*invJ(3,3)*Zeta2;
                B(5,25) = 0;
                B(0,26) = 0;
                B(1,26) = 0;
                B(2,26) = 4*invJ(1,3)*Zeta4+4*invJ(3,3)*Zeta2;
                B(3,26) = 0;
                B(4,26) = 4*invJ(1,2)*Zeta4+4*invJ(3,2)*Zeta2;
                B(5,26) = 4*invJ(1,1)*Zeta4+4*invJ(3,1)*Zeta2;
                B(0,27) = 4*invJ(2,1)*Zeta4+4*invJ(3,1)*Zeta3;
                B(1,27) = 0;
                B(2,27) = 0;
                B(3,27) = 4*invJ(2,2)*Zeta4+4*invJ(3,2)*Zeta3;
                B(4,27) = 0;
                B(5,27) = 4*invJ(2,3)*Zeta4+4*invJ(3,3)*Zeta3;
                B(0,28) = 0;
                B(1,28) = 4*invJ(2,2)*Zeta4+4*invJ(3,2)*Zeta3;
                B(2,28) = 0;
                B(3,28) = 4*invJ(2,1)*Zeta4+4*invJ(3,1)*Zeta3;
                B(4,28) = 4*invJ(2,3)*Zeta4+4*invJ(3,3)*Zeta3;
                B(5,28) = 0;
                B(0,29) = 0;
                B(1,29) = 0;
                B(2,29) = 4*invJ(2,3)*Zeta4+4*invJ(3,3)*Zeta3;
                B(3,29) = 0;
                B(4,29) = 4*invJ(2,2)*Zeta4+4*invJ(3,2)*Zeta3;
                B(5,29) = 4*invJ(2,1)*Zeta4+4*invJ(3,1)*Zeta3;
        }

        virtual Matrix InterpGaussToNodal(MatrixType eMatType, const Matrix& taug) const {
                const IntegrationRule& oIntegRule = GetIntegrationRule(eMatType);
                const octave_idx_type iNumGauss = oIntegRule.iGetNumEvalPoints();
                const octave_idx_type iNumDir = oIntegRule.iGetNumDirections();
                const octave_idx_type iNumNodesCorner = 4;
                ColumnVector rv(iNumDir);

                Matrix H(iNumGauss, iNumNodesCorner);

                for (octave_idx_type i = 0; i < iNumGauss; ++i) {
                        for (octave_idx_type j = 0; j < iNumDir; ++j) {
                                rv(j) = oIntegRule.dGetPosition(i, j);
                        }

                        StressInterpMatrix(rv, H, i);
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
                                taun(idx.imid, j) = 0.5 * (taun(idx.ico1, j) + taun(idx.ico2, j));
                        }
                }

                return taun;
        }

private:
        void StressInterpMatrix(const ColumnVector& rv, Matrix& Hs, octave_idx_type irow) const {
                FEM_ASSERT(rv.numel() == 4);
                FEM_ASSERT(Hs.columns() == 4);
                FEM_ASSERT(irow >= 0);
                FEM_ASSERT(irow < Hs.rows());

                const double Zeta1 = rv(0);
                const double Zeta2 = rv(1);
                const double Zeta3 = rv(2);
                const double Zeta4 = rv(3);

                Hs(irow, 0) = Zeta1;
                Hs(irow, 1) = Zeta2;
                Hs(irow, 2) = Zeta3;
                Hs(irow, 3) = Zeta4;
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
                return rv(0) + rv(1) + rv(2) - 1;
        }

        static void GetElemLimits(ColumnVector& rmin, ColumnVector& rmax) {
                FEM_ASSERT(rmin.rows() == iGetNumDirections());
                FEM_ASSERT(rmax.rows() == rmin.rows());

                for (octave_idx_type i = 0; i < rmin.rows(); ++i) {
                        rmin(i) = 0.;
                }

                for (octave_idx_type i = 0; i < rmax.rows(); ++i) {
                        rmax(i) = 1.;
                }
        }

        static void ScalarInterpMatrix(const ColumnVector& rv, RowVector& HA) {
                FEM_ASSERT(rv.rows() == 3);
                FEM_ASSERT(HA.columns() == 6);

                const double Zeta1 = rv(0);
                const double Zeta2 = rv(1);
                const double Zeta3 = rv(2);

                HA(0) = Zeta1*(2*Zeta1-1);
                HA(1) = Zeta2*(2*Zeta2-1);
                HA(2) = Zeta3*(2*Zeta3-1);
                HA(3) = 4*Zeta1*Zeta2;
                HA(4) = 4*Zeta2*Zeta3;
                HA(5) = 4*Zeta1*Zeta3;
        }

        static void VectorInterpMatrix(const ColumnVector& rv, Matrix& Hf) {
                FEM_ASSERT(rv.rows() == 3);
                FEM_ASSERT(Hf.rows() == 3);
                FEM_ASSERT(Hf.columns() == 18);

                const double Zeta1 = rv(0);
                const double Zeta2 = rv(1);
                const double Zeta3 = rv(2);

                Hf(0,0) = Zeta1*(2*Zeta1-1);
                Hf(1,0) = 0;
                Hf(2,0) = 0;
                Hf(0,1) = 0;
                Hf(1,1) = Zeta1*(2*Zeta1-1);
                Hf(2,1) = 0;
                Hf(0,2) = 0;
                Hf(1,2) = 0;
                Hf(2,2) = Zeta1*(2*Zeta1-1);
                Hf(0,3) = Zeta2*(2*Zeta2-1);
                Hf(1,3) = 0;
                Hf(2,3) = 0;
                Hf(0,4) = 0;
                Hf(1,4) = Zeta2*(2*Zeta2-1);
                Hf(2,4) = 0;
                Hf(0,5) = 0;
                Hf(1,5) = 0;
                Hf(2,5) = Zeta2*(2*Zeta2-1);
                Hf(0,6) = Zeta3*(2*Zeta3-1);
                Hf(1,6) = 0;
                Hf(2,6) = 0;
                Hf(0,7) = 0;
                Hf(1,7) = Zeta3*(2*Zeta3-1);
                Hf(2,7) = 0;
                Hf(0,8) = 0;
                Hf(1,8) = 0;
                Hf(2,8) = Zeta3*(2*Zeta3-1);
                Hf(0,9) = 4*Zeta1*Zeta2;
                Hf(1,9) = 0;
                Hf(2,9) = 0;
                Hf(0,10) = 0;
                Hf(1,10) = 4*Zeta1*Zeta2;
                Hf(2,10) = 0;
                Hf(0,11) = 0;
                Hf(1,11) = 0;
                Hf(2,11) = 4*Zeta1*Zeta2;
                Hf(0,12) = 4*Zeta2*Zeta3;
                Hf(1,12) = 0;
                Hf(2,12) = 0;
                Hf(0,13) = 0;
                Hf(1,13) = 4*Zeta2*Zeta3;
                Hf(2,13) = 0;
                Hf(0,14) = 0;
                Hf(1,14) = 0;
                Hf(2,14) = 4*Zeta2*Zeta3;
                Hf(0,15) = 4*Zeta1*Zeta3;
                Hf(1,15) = 0;
                Hf(2,15) = 0;
                Hf(0,16) = 0;
                Hf(1,16) = 4*Zeta1*Zeta3;
                Hf(2,16) = 0;
                Hf(0,17) = 0;
                Hf(1,17) = 0;
                Hf(2,17) = 4*Zeta1*Zeta3;
        }

        static void VectorInterpMatrixDerR(const ColumnVector& rv, Matrix& dHf_dr) {
                FEM_ASSERT(rv.rows() == 3);
                FEM_ASSERT(dHf_dr.rows() == 3);
                FEM_ASSERT(dHf_dr.columns() == 18);

                const double Zeta1 = rv(0);
                const double Zeta2 = rv(1);
                const double Zeta3 = rv(2);

                dHf_dr(0,0) = 4*Zeta1-1;
                dHf_dr(1,0) = 0;
                dHf_dr(2,0) = 0;
                dHf_dr(0,1) = 0;
                dHf_dr(1,1) = 4*Zeta1-1;
                dHf_dr(2,1) = 0;
                dHf_dr(0,2) = 0;
                dHf_dr(1,2) = 0;
                dHf_dr(2,2) = 4*Zeta1-1;
                dHf_dr(0,3) = 0;
                dHf_dr(1,3) = 0;
                dHf_dr(2,3) = 0;
                dHf_dr(0,4) = 0;
                dHf_dr(1,4) = 0;
                dHf_dr(2,4) = 0;
                dHf_dr(0,5) = 0;
                dHf_dr(1,5) = 0;
                dHf_dr(2,5) = 0;
                dHf_dr(0,6) = 1-4*Zeta3;
                dHf_dr(1,6) = 0;
                dHf_dr(2,6) = 0;
                dHf_dr(0,7) = 0;
                dHf_dr(1,7) = 1-4*Zeta3;
                dHf_dr(2,7) = 0;
                dHf_dr(0,8) = 0;
                dHf_dr(1,8) = 0;
                dHf_dr(2,8) = 1-4*Zeta3;
                dHf_dr(0,9) = 4*Zeta2;
                dHf_dr(1,9) = 0;
                dHf_dr(2,9) = 0;
                dHf_dr(0,10) = 0;
                dHf_dr(1,10) = 4*Zeta2;
                dHf_dr(2,10) = 0;
                dHf_dr(0,11) = 0;
                dHf_dr(1,11) = 0;
                dHf_dr(2,11) = 4*Zeta2;
                dHf_dr(0,12) = -4*Zeta2;
                dHf_dr(1,12) = 0;
                dHf_dr(2,12) = 0;
                dHf_dr(0,13) = 0;
                dHf_dr(1,13) = -4*Zeta2;
                dHf_dr(2,13) = 0;
                dHf_dr(0,14) = 0;
                dHf_dr(1,14) = 0;
                dHf_dr(2,14) = -4*Zeta2;
                dHf_dr(0,15) = 4*Zeta3-4*Zeta1;
                dHf_dr(1,15) = 0;
                dHf_dr(2,15) = 0;
                dHf_dr(0,16) = 0;
                dHf_dr(1,16) = 4*Zeta3-4*Zeta1;
                dHf_dr(2,16) = 0;
                dHf_dr(0,17) = 0;
                dHf_dr(1,17) = 0;
                dHf_dr(2,17) = 4*Zeta3-4*Zeta1;
        }

        static void VectorInterpMatrixDerS(const ColumnVector& rv, Matrix& dHf_ds) {
                FEM_ASSERT(rv.rows() == 3);
                FEM_ASSERT(dHf_ds.rows() == 3);
                FEM_ASSERT(dHf_ds.columns() == 18);

                const double Zeta1 = rv(0);
                const double Zeta2 = rv(1);
                const double Zeta3 = rv(2);

                dHf_ds(0,0) = 0;
                dHf_ds(1,0) = 0;
                dHf_ds(2,0) = 0;
                dHf_ds(0,1) = 0;
                dHf_ds(1,1) = 0;
                dHf_ds(2,1) = 0;
                dHf_ds(0,2) = 0;
                dHf_ds(1,2) = 0;
                dHf_ds(2,2) = 0;
                dHf_ds(0,3) = 4*Zeta2-1;
                dHf_ds(1,3) = 0;
                dHf_ds(2,3) = 0;
                dHf_ds(0,4) = 0;
                dHf_ds(1,4) = 4*Zeta2-1;
                dHf_ds(2,4) = 0;
                dHf_ds(0,5) = 0;
                dHf_ds(1,5) = 0;
                dHf_ds(2,5) = 4*Zeta2-1;
                dHf_ds(0,6) = 1-4*Zeta3;
                dHf_ds(1,6) = 0;
                dHf_ds(2,6) = 0;
                dHf_ds(0,7) = 0;
                dHf_ds(1,7) = 1-4*Zeta3;
                dHf_ds(2,7) = 0;
                dHf_ds(0,8) = 0;
                dHf_ds(1,8) = 0;
                dHf_ds(2,8) = 1-4*Zeta3;
                dHf_ds(0,9) = 4*Zeta1;
                dHf_ds(1,9) = 0;
                dHf_ds(2,9) = 0;
                dHf_ds(0,10) = 0;
                dHf_ds(1,10) = 4*Zeta1;
                dHf_ds(2,10) = 0;
                dHf_ds(0,11) = 0;
                dHf_ds(1,11) = 0;
                dHf_ds(2,11) = 4*Zeta1;
                dHf_ds(0,12) = 4*Zeta3-4*Zeta2;
                dHf_ds(1,12) = 0;
                dHf_ds(2,12) = 0;
                dHf_ds(0,13) = 0;
                dHf_ds(1,13) = 4*Zeta3-4*Zeta2;
                dHf_ds(2,13) = 0;
                dHf_ds(0,14) = 0;
                dHf_ds(1,14) = 0;
                dHf_ds(2,14) = 4*Zeta3-4*Zeta2;
                dHf_ds(0,15) = -4*Zeta1;
                dHf_ds(1,15) = 0;
                dHf_ds(2,15) = 0;
                dHf_ds(0,16) = 0;
                dHf_ds(1,16) = -4*Zeta1;
                dHf_ds(2,16) = 0;
                dHf_ds(0,17) = 0;
                dHf_ds(1,17) = 0;
                dHf_ds(2,17) = -4*Zeta1;
        }

        static const IntegrationRule& GetIntegrationRule(Element::MatrixType eMatType) {
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
                case Element::VEC_LOAD_CONSISTENT: {
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
                        rmin(i) = -1.;
                }

                for (octave_idx_type i = 0; i < rmax.rows(); ++i) {
                        rmax(i) = 1.;
                }
        }

        static void ScalarInterpMatrix(const ColumnVector& rv, RowVector& HA) {
                FEM_ASSERT(rv.rows() == 2);

                const double r = rv(0);
                const double s = rv(1);

                HA(0) = ((r+1)*(s+1))/4.0;
                HA(1) = ((1-r)*(s+1))/4.0;
                HA(2) = ((1-r)*(1-s))/4.0;
                HA(3) = ((r+1)*(1-s))/4.0;
        }

        static void VectorInterpMatrix(const ColumnVector& rv, Matrix& Hf) {
                FEM_ASSERT(rv.rows() == 2);

                const double r = rv(0);
                const double s = rv(1);

                Hf(0,0) = ((r+1)*(s+1))/4.0;
                Hf(1,0) = 0;
                Hf(2,0) = 0;
                Hf(0,1) = 0;
                Hf(1,1) = ((r+1)*(s+1))/4.0;
                Hf(2,1) = 0;
                Hf(0,2) = 0;
                Hf(1,2) = 0;
                Hf(2,2) = ((r+1)*(s+1))/4.0;
                Hf(0,3) = ((1-r)*(s+1))/4.0;
                Hf(1,3) = 0;
                Hf(2,3) = 0;
                Hf(0,4) = 0;
                Hf(1,4) = ((1-r)*(s+1))/4.0;
                Hf(2,4) = 0;
                Hf(0,5) = 0;
                Hf(1,5) = 0;
                Hf(2,5) = ((1-r)*(s+1))/4.0;
                Hf(0,6) = ((1-r)*(1-s))/4.0;
                Hf(1,6) = 0;
                Hf(2,6) = 0;
                Hf(0,7) = 0;
                Hf(1,7) = ((1-r)*(1-s))/4.0;
                Hf(2,7) = 0;
                Hf(0,8) = 0;
                Hf(1,8) = 0;
                Hf(2,8) = ((1-r)*(1-s))/4.0;
                Hf(0,9) = ((r+1)*(1-s))/4.0;
                Hf(1,9) = 0;
                Hf(2,9) = 0;
                Hf(0,10) = 0;
                Hf(1,10) = ((r+1)*(1-s))/4.0;
                Hf(2,10) = 0;
                Hf(0,11) = 0;
                Hf(1,11) = 0;
                Hf(2,11) = ((r+1)*(1-s))/4.0;

        }

        static void VectorInterpMatrixDerR(const ColumnVector& rv, Matrix& dHf_dr) {
                FEM_ASSERT(rv.rows() == 2);

                // const double r = rv(0);
                const double s = rv(1);

                dHf_dr(0,0) = (s+1)/4.0;
                dHf_dr(1,0) = 0;
                dHf_dr(2,0) = 0;
                dHf_dr(0,1) = 0;
                dHf_dr(1,1) = (s+1)/4.0;
                dHf_dr(2,1) = 0;
                dHf_dr(0,2) = 0;
                dHf_dr(1,2) = 0;
                dHf_dr(2,2) = (s+1)/4.0;
                dHf_dr(0,3) = -(s+1)/4.0;
                dHf_dr(1,3) = 0;
                dHf_dr(2,3) = 0;
                dHf_dr(0,4) = 0;
                dHf_dr(1,4) = -(s+1)/4.0;
                dHf_dr(2,4) = 0;
                dHf_dr(0,5) = 0;
                dHf_dr(1,5) = 0;
                dHf_dr(2,5) = -(s+1)/4.0;
                dHf_dr(0,6) = -(1-s)/4.0;
                dHf_dr(1,6) = 0;
                dHf_dr(2,6) = 0;
                dHf_dr(0,7) = 0;
                dHf_dr(1,7) = -(1-s)/4.0;
                dHf_dr(2,7) = 0;
                dHf_dr(0,8) = 0;
                dHf_dr(1,8) = 0;
                dHf_dr(2,8) = -(1-s)/4.0;
                dHf_dr(0,9) = (1-s)/4.0;
                dHf_dr(1,9) = 0;
                dHf_dr(2,9) = 0;
                dHf_dr(0,10) = 0;
                dHf_dr(1,10) = (1-s)/4.0;
                dHf_dr(2,10) = 0;
                dHf_dr(0,11) = 0;
                dHf_dr(1,11) = 0;
                dHf_dr(2,11) = (1-s)/4.0;
        }

        static void VectorInterpMatrixDerS(const ColumnVector& rv, Matrix& dHf_ds) {
                FEM_ASSERT(rv.rows() == 2);

                const double r = rv(0);
                // const double s = rv(1);

                dHf_ds(0,0) = (r+1)/4.0;
                dHf_ds(1,0) = 0;
                dHf_ds(2,0) = 0;
                dHf_ds(0,1) = 0;
                dHf_ds(1,1) = (r+1)/4.0;
                dHf_ds(2,1) = 0;
                dHf_ds(0,2) = 0;
                dHf_ds(1,2) = 0;
                dHf_ds(2,2) = (r+1)/4.0;
                dHf_ds(0,3) = (1-r)/4.0;
                dHf_ds(1,3) = 0;
                dHf_ds(2,3) = 0;
                dHf_ds(0,4) = 0;
                dHf_ds(1,4) = (1-r)/4.0;
                dHf_ds(2,4) = 0;
                dHf_ds(0,5) = 0;
                dHf_ds(1,5) = 0;
                dHf_ds(2,5) = (1-r)/4.0;
                dHf_ds(0,6) = -(1-r)/4.0;
                dHf_ds(1,6) = 0;
                dHf_ds(2,6) = 0;
                dHf_ds(0,7) = 0;
                dHf_ds(1,7) = -(1-r)/4.0;
                dHf_ds(2,7) = 0;
                dHf_ds(0,8) = 0;
                dHf_ds(1,8) = 0;
                dHf_ds(2,8) = -(1-r)/4.0;
                dHf_ds(0,9) = -(r+1)/4.0;
                dHf_ds(1,9) = 0;
                dHf_ds(2,9) = 0;
                dHf_ds(0,10) = 0;
                dHf_ds(1,10) = -(r+1)/4.0;
                dHf_ds(2,10) = 0;
                dHf_ds(0,11) = 0;
                dHf_ds(1,11) = 0;
                dHf_ds(2,11) = -(r+1)/4.0;
        }

        static const IntegrationRule& GetIntegrationRule(Element::MatrixType eMatType) {
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

class PressureLoad: public Element {
public:
        PressureLoad(octave_idx_type id, const Matrix& X, const Material* material, const int32NDArray& nodes, const Matrix& p, octave_idx_type colidx)
                :Element(id, X, material, nodes), p(p), colidx(colidx) {

                FEM_ASSERT(X.rows() == 3);
                FEM_ASSERT(X.columns() == p.columns());
                FEM_ASSERT(X.columns() == nodes.numel());
        }

        virtual void Assemble(MatrixAss& mat, MeshInfo& info, const DofMap& dof, MatrixType eMatType) const {
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
                                rv(j) = oIntegRule.dGetPosition(i, j);
                        }

                        DisplacementInterpMatrix(rv, Hf);
                        DisplacementInterpMatrixDerR(rv, dHf_dr);
                        DisplacementInterpMatrixDerS(rv, dHf_ds);
                        PressureInterpMatrix(rv, HA);

                        SurfaceNormalVector(dHf_dr, n1);
                        SurfaceNormalVector(dHf_ds, n2);

                        n_detJA(0) = n1(1) * n2(2) - n1(2) * n2(1);
                        n_detJA(1) = n1(2) * n2(0) - n1(0) * n2(2);
                        n_detJA(2) = n1(0) * n2(1) - n1(1) * n2(0);

                        double detJA_2 = 0.;

                        for (octave_idx_type l = 0; l < 3; ++l) {
                                detJA_2 += n_detJA(l) * n_detJA(l);
                        }

                        for (octave_idx_type l = 0; l < iNumDof; ++l) {
                                double HfT_nl_detJA = 0.;

                                for (octave_idx_type m = 0; m < 3; ++m) {
                                        HfT_nl_detJA -= Hf(m, l) * n_detJA(m);
                                }

                                HfT_n_dA(l) = HfT_nl_detJA * alpha;
                        }

                        for (octave_idx_type l = 0; l < iNumLoads; ++l) {
                                double HA_pl = 0.;

                                for (octave_idx_type m = 0; m < iNumNodes; ++m) {
                                        HA_pl += HA(m) * p(l, m);
                                }

                                HA_p(l) = HA_pl;
                        }

                        for (octave_idx_type l = 0; l < iNumLoads; ++l) {
                                for (octave_idx_type m = 0; m < iNumDof; ++m) {
                                        fA(m, l) += HfT_n_dA(m) * HA_p(l);
                                }
                        }
                }

                for (octave_idx_type j = 0; j < iNumLoads; ++j) {
                        for (octave_idx_type i = 0; i < iNumDof; ++i) {
                                mat.Insert(fA(i, j), dofidx(i), colidx + j);
                        }
                }
        }

        virtual octave_idx_type iGetWorkSpaceSize(MatrixType eMatType) const {
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

protected:
        void SurfaceNormalVector(const Matrix& dHf, ColumnVector& n) const {
                for (octave_idx_type i = 0; i < 3; ++i) {
                        double ni = 0.;

                        for (octave_idx_type j = 0; j < nodes.numel(); ++j) {
                                for (octave_idx_type k = 0; k < 3; ++k) {
                                        ni += dHf(i, j * 3 + k) * X(k, j);
                                }
                        }

                        n(i) = ni;
                }
        }
        virtual const IntegrationRule& GetIntegrationRule(MatrixType eMatType) const=0;
        virtual void PressureInterpMatrix(const ColumnVector& r, RowVector& HA) const=0;
        virtual void DisplacementInterpMatrix(const ColumnVector& r, Matrix& Hf) const=0;
        virtual void DisplacementInterpMatrixDerR(const ColumnVector& r, Matrix& dHf_dr) const=0;
        virtual void DisplacementInterpMatrixDerS(const ColumnVector& r, Matrix& dHf_ds) const=0;

private:
        Matrix p;
        octave_idx_type colidx;
};

template <typename SHAPE_FUNC>
class PressureLoadImp: public PressureLoad {
public:
        PressureLoadImp(octave_idx_type id, const Matrix& X, const Material* material, const int32NDArray& nodes, const Matrix& p, octave_idx_type colidx)
                :PressureLoad(id, X, material, nodes, p, colidx) {
                FEM_ASSERT(nodes.numel() == SHAPE_FUNC::iGetNumNodes());
        }

protected:
        void PressureInterpMatrix(const ColumnVector& rv, RowVector& HA) const {
                SHAPE_FUNC::ScalarInterpMatrix(rv, HA);
        }

        void DisplacementInterpMatrix(const ColumnVector& rv, Matrix& Hf) const {
                SHAPE_FUNC::VectorInterpMatrix(rv, Hf);
        }

        void DisplacementInterpMatrixDerR(const ColumnVector& rv, Matrix& dHf_dr) const {
                SHAPE_FUNC::VectorInterpMatrixDerR(rv, dHf_dr);
        }

        void DisplacementInterpMatrixDerS(const ColumnVector& rv, Matrix& dHf_ds) const {
                SHAPE_FUNC::VectorInterpMatrixDerS(rv, dHf_ds);
        }

        virtual const IntegrationRule& GetIntegrationRule(MatrixType eMatType) const {
                return SHAPE_FUNC::GetIntegrationRule(eMatType);
        }
};

typedef PressureLoadImp<ShapeIso4> PressureLoadIso4;
typedef PressureLoadImp<ShapeTria6> PressureLoadTria6;

class StructForce: public Element {
public:
        StructForce(octave_idx_type id, const Matrix& X, const Material* material, const int32NDArray& nodes, octave_idx_type colidx, const Matrix& loads)
                :Element(id, X, material, nodes),
                 loads(loads),
                 colidx(colidx) {

                FEM_ASSERT(X.rows() == 3);
                FEM_ASSERT(loads.columns() == 3 || loads.columns() == 6);
                FEM_ASSERT(loads.rows() == nodes.rows());
        }

        virtual void Assemble(MatrixAss& mat, MeshInfo& info, const DofMap& dof, MatrixType eMatType) const {
                switch (eMatType) {
                case VEC_LOAD_CONSISTENT:
                case VEC_LOAD_LUMPED:
                        break;
                default:
                        return;
                }

                for (octave_idx_type j = 0; j < loads.columns(); ++j) {
                        for (octave_idx_type i = 0; i < loads.rows(); ++i) {
                                const octave_idx_type inode = nodes(i).value() - 1;

                                mat.Insert(loads(i, j), dof.GetNodeDofIndex(inode, j), colidx);
                        }
                }
        }

        virtual octave_idx_type iGetWorkSpaceSize(MatrixType eMatType) const {
                switch (eMatType) {
                case VEC_LOAD_CONSISTENT:
                case VEC_LOAD_LUMPED:
                        return loads.rows() * loads.columns();
                default:
                        return 0;
                }
        }

private:
        Matrix loads;
        octave_idx_type colidx;
};

class ElementTypes {
public:
        enum TypeId {
                ELEM_ISO8 = 0,
                ELEM_TET10,
                ELEM_RBE3,
                ELEM_JOINT,
                ELEM_SFNCON4,
                ELEM_SFNCON6,
                ELEM_PRESSURE_ISO4,
                ELEM_PRESSURE_TRIA6,
                ELEM_STRUCT_FORCE,
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
        {ElementTypes::ELEM_ISO8,           "iso8",     8,  8, DofMap::ELEM_NODOF},
        {ElementTypes::ELEM_TET10,          "tet10",   10, 10, DofMap::ELEM_NODOF},
        {ElementTypes::ELEM_RBE3,           "rbe3",     2, -1, DofMap::ELEM_RBE3},
        {ElementTypes::ELEM_JOINT,          "joints",   1, -1, DofMap::ELEM_JOINT},
        {ElementTypes::ELEM_SFNCON4,        "sfncon4",  1, -1, DofMap::ELEM_JOINT},
        {ElementTypes::ELEM_SFNCON6,        "sfncon6",  1, -1, DofMap::ELEM_JOINT},
        {ElementTypes::ELEM_PRESSURE_ISO4,  "iso4",     4,  4, DofMap::ELEM_NODOF},
        {ElementTypes::ELEM_PRESSURE_TRIA6, "tria6",    6,  6, DofMap::ELEM_NODOF},
        {ElementTypes::ELEM_STRUCT_FORCE,   "force",    1, -1, DofMap::ELEM_NODOF}
};

class ElementBlockBase {
public:
        explicit ElementBlockBase(ElementTypes::TypeId eltype)
                :eltype(eltype) {
        }

        virtual ~ElementBlockBase() {
        }

        virtual octave_idx_type iGetWorkSpaceSize(Element::MatrixType eMatType) const=0;
        virtual void Assemble(MatrixAss& oMatAss, MeshInfo& info, const DofMap& oDof, Element::MatrixType eMatType) const=0;
        virtual void PostProcElem(NDArray& mat, Element::MatrixType eMatType, const NDArray& U) const=0;
        virtual double dGetMass() const=0;
        virtual bool bNeedMatrixInfo(Element::MatrixType eMatType) const=0;

        template <typename ElementType, typename... Args>
        inline void Insert(octave_idx_type id, const Matrix& X, const Material* material, const int32NDArray& nodes, const Args&... args);

        ElementTypes::TypeId GetElementType() const { return eltype; }
        virtual void Extract(octave_idx_type& idx, octave_map& sElem) const=0;
        virtual octave_idx_type iGetNumElem() const=0;

private:
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
                        for (octave_idx_type j = 0; j < elements.columns(); ++j) {
                                nodes_e(j) = elements(i, j);
                                for (octave_idx_type k = 0; k < X_e.rows(); ++k) {
                                        X_e(k, j) = nodes(nodes_e(j).value() - 1, k);
                                }
                        }

                        const Material* material = nullptr; // Some elements like RBE3 do not need a material

                        if (materials(i).value() > 0) {
                                material = &rgMaterials[materials(i).value() - 1];
                        }

                        rgElements.emplace_back(i + 1, X_e, material, nodes_e, args...);
                }
        }

        template <typename... Args>
        void Insert(octave_idx_type id, const Matrix& X, const Material* material, const int32NDArray& nodes, const Args&... args) {
                rgElements.emplace_back(id, X, material, nodes, args...);
        }

        octave_idx_type iGetWorkSpaceSize(Element::MatrixType eMatType) const {
                octave_idx_type iWorkSpace = 0;

                for (auto i = rgElements.begin(); i != rgElements.end(); ++i) {
                        iWorkSpace += i->iGetWorkSpaceSize(eMatType);
                }

                return iWorkSpace;
        }

        void Assemble(MatrixAss& oMatAss, MeshInfo& oMeshInfo, const DofMap& oDof, Element::MatrixType eMatType) const {
                for (auto i = rgElements.begin(); i != rgElements.end(); ++i) {
                        i->Assemble(oMatAss, oMeshInfo, oDof, eMatType);

                        OCTAVE_QUIT;
                }
        }

        void PostProcElem(NDArray& mat, Element::MatrixType eMatType, const NDArray& U) const {
                for (auto i = rgElements.begin(); i != rgElements.end(); ++i) {
                        i->PostProcElem(mat, eMatType, U);

                        OCTAVE_QUIT;
                }
        }

        double dGetMass() const {
                double dm = 0.;

                for (auto i = rgElements.begin(); i != rgElements.end(); ++i) {
                        dm += i->dGetMass();

                        OCTAVE_QUIT;
                }

                return dm;
        }

        virtual bool bNeedMatrixInfo(Element::MatrixType eMatType) const {
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
                        oElem.Extract(idx, sElem);
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

                                        const octave_scalar_map pressure = cell_pressure(i).scalar_map_value();

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
                                                throw std::runtime_error("pressure.elements number of columns does not match");
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
                                                        for (octave_idx_type l = 0; l < X.columns(); ++l) {
                                                                for (octave_idx_type m = 0; m < X.rows(); ++m) {
                                                                        octave_idx_type inode = elements(k, l).value() - 1;

                                                                        if (inode < 0 || inode >= nodes.rows()) {
                                                                                throw std::runtime_error("node index out of range in pressure.elements");
                                                                        }

                                                                        X(m, l) = nodes(inode, m);
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

namespace shape_func_util {
        template <typename T>
        struct SelectElemPerSlaveNode {
        };

        template <>
        struct SelectElemPerSlaveNode<ShapeIso4> {
                static constexpr octave_idx_type iNumElemPerSlaveNode = 9;
        };

        template <>
        struct SelectElemPerSlaveNode<ShapeTria6> {
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
                        throw std::runtime_error("invalid value for elements.sfncon{4|6}.constraint");
                }
        }

        static ConstraintType GetConstraintType(const Cell& ov, octave_idx_type j) {
                if (!ov.numel()) {
                        return CT_FIXED;
                }

                const int constr = ov(j).int_value();

                if (error_state) {
                        throw std::runtime_error("elements.sfncon{4|6}.constraint must be a scalar value");
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
const char SurfToNodeConstrBase::szErrCompileWithNlopt[] = "__mboct_fem_pkg__ must be compiled with nlopt in order to use element types sfncon{4|6}";
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
                                Xs(l) = X(nidxslave(k).value() - 1, l);
                        }

                        for (octave_idx_type i = 0; i < nidxmaster.rows(); ++i) {
                                double dXmin = std::numeric_limits<double>::max();

                                for (octave_idx_type j = 0; j < nidxmaster.columns(); ++j) {
                                        double dX = 0;

                                        for (octave_idx_type l = 0; l < Xs.rows(); ++l) {
                                                dX += std::pow(X(nidxmaster(i, j).value() - 1, l) - Xs(l), 2);
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
                                Xs(j) = X(nidxslave(i).value() - 1, j);
                        }

                        double fopt = std::numeric_limits<double>::max();
                        nlopt_result rcopt = NLOPT_FAILURE;
                        octave_idx_type lopt = -1;
                        rvopt.fill(0.);

                        for (octave_idx_type l = 0; l < eidxmaster[i].size(); ++l) {
                                for (octave_idx_type j = 0; j < nidxmaster.columns(); ++j) {
                                        for (octave_idx_type k = 0; k < iNumDofNode; ++k) {
                                                Xm(j * iNumDofNode + k) = X(nidxmaster(eidxmaster[i][l].eidx, j).value() - 1, k);
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
                                                os << nidxmaster(eidx, k).value() << ' ';
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
                                        os << X(nidxslave(i).value() - 1, k) << ' ';
                                }

                                os << std::endl;

                                os << "Xi=" << Xiopt.transpose() << std::endl;
                                os << "Xm=" << std::endl;

                                if (eidx >= 0) {
                                        for (octave_idx_type k = 0; k < iNumNodesElem; ++k) {
                                                for (octave_idx_type j = 0; j < iNumDofNode; ++j) {
                                                        os << X(nidxmaster(eidx, k).value() - 1, j) << " ";
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

                        C.fill(0.);

                        for (octave_idx_type j = 0; j < 6; ++j) {
                                for (octave_idx_type k = 0; k < C.rows(); ++k) {
                                        C(k, j) = (k == j);
                                }
                        }

                        for (octave_idx_type j = 0; j < nidxmaster.columns(); ++j) {
                                for (octave_idx_type k = 0; k < iNumDofNode; ++k) {
                                        for (octave_idx_type l = 0; l < C.rows(); ++l) {
                                                C(l, (j + 1) * 6 + k) = -Hf(l, j * iNumDofNode + k);
                                        }
                                }
                        }

                        enodes(0) = nidxslave(i);

                        for (octave_idx_type j = 0; j < nidxmaster.columns(); ++j) {
                                enodes(j + 1) = nidxmaster(eidxmaster[i][lopt].eidx, j).value();
                        }

                        for (octave_idx_type j = 0; j < X.columns(); ++j) {
                                Xe(j, 0) = X(nidxslave(i).value() - 1, j);
                        }

                        for (octave_idx_type j = 0; j < nidxmaster.columns(); ++j) {
                                for (octave_idx_type k = 0; k < X.columns(); ++k) {
                                        Xe(k, j + 1) = X(nidxmaster(eidxmaster[i][lopt].eidx, j).value() - 1, k);
                                }
                        }

                        switch (eType) {
                        case CT_SLIDING: {
                                SHAPE_FUNC::VectorInterpMatrixDerR(rvopt, dHf_dr);
                                SHAPE_FUNC::VectorInterpMatrixDerS(rvopt, dHf_ds);
                                SurfaceNormalVector(Xe, dHf_dr, n1);
                                SurfaceNormalVector(Xe, dHf_ds, n2);

                                n(0) = n1(1) * n2(2) - n1(2) * n2(1);
                                n(1) = n1(2) * n2(0) - n1(0) * n2(2);
                                n(2) = n1(0) * n2(1) - n1(1) * n2(0);

                                double norm_n = 0.;

                                for (octave_idx_type l = 0; l < 3; ++l) {
                                        norm_n += n(l) * n(l);
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
                                        n(l) /= norm_n;
                                }

                                nC.fill(0.);

                                for (octave_idx_type j = 0; j < C.columns(); ++j) {
                                        for (octave_idx_type k = 0; k < C.rows(); ++k) {
                                                nC(j) += n(k) * C(k, j);
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
                                        ni += dHf(i, j * 3 + k) * X(k, j + 1);
                                }
                        }

                        n(i) = ni;
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
                        dX(i, 0) = std::numeric_limits<double>::max();
                        dX(i, 1) = -dX(i, 0);

                        for (octave_idx_type j = 0; j < iNumNodesElem; ++j) {
                                dX(i, 0) = std::min(dX(i, 0), Xm(j * iNumDofNode + i));
                                dX(i, 1) = std::max(dX(i, 1), Xm(j * iNumDofNode + i));
                        }
                }

                double dTolF = 0;

                for (octave_idx_type i = 0; i < iNumDofNode; ++i) {
                        dTolF = std::max(dTolF, dTolX * (dX(i, 1) - dX(i, 0)));
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

                FEM_ASSERT(n == pFuncData->rv.rows());

                for (octave_idx_type i = 0; i < pFuncData->rv.rows(); ++i) {
                        pFuncData->rv(i) = x[i];
                }

                return pFuncData->oSNCO.Objective(pFuncData->rv, pFuncData->f, pFuncData->Hf);
        }

        static double EqualityConstr(unsigned n, const double x[], double gradient[], void* pData) {
                auto pFuncData = static_cast<FuncData*>(pData);

                FEM_ASSERT(n == pFuncData->rv.rows());

                for (octave_idx_type i = 0; i < pFuncData->rv.rows(); ++i) {
                        pFuncData->rv(i) = x[i];
                }

                return pFuncData->oSNCO.EqualityConstr(pFuncData->rv);
        }

        double Objective(const ColumnVector& rv, ColumnVector& f, Matrix& Hf) const {
                FEM_ASSERT(f.rows() == iNumDofNode);
                FEM_ASSERT(f.rows() == Xs.rows());

                Position(rv, f, Hf);

                for (octave_idx_type i = 0; i < f.rows(); ++i) {
                        f(i) -= Xs(i);
                }

                double ftot = 0.;

                for (octave_idx_type i = 0; i < f.rows(); ++i) {
                        ftot += f(i) * f(i);
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
                        Xi(i) = 0;
                }

                for (octave_idx_type j = 0; j < Hf.columns(); ++j) {
                        for (octave_idx_type i = 0; i < Hf.rows(); ++i) {
                                Xi(i) += Hf(i, j) * Xm(j);
                        }
                }
        }

private:
        const ColumnVector Xm, Xs;

        struct FuncData {
                FuncData(const ColumnVector& Xm, const ColumnVector& Xs)
                        :oSNCO(Xm, Xs),
                         opt(nlopt_create(NLOPT_LN_COBYLA, iNumDir)),
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
                throw std::runtime_error("elements.sfncon{4|6} must be a struct array");
        }

        const auto iter_nidxmaster = s_elem.seek("master");

        if (iter_nidxmaster == s_elem.end()) {
                throw std::runtime_error("elements.sfncon{4|6}.master not defined");
        }

        const Cell ov_nidxmaster = s_elem.contents(iter_nidxmaster);

        const auto iter_nidxslave = s_elem.seek("slave");

        if (iter_nidxslave == s_elem.end()) {
                throw std::runtime_error("elements.sfncon{4|6}.slave not defined");
        }

        const Cell ov_nidxslave = s_elem.contents(iter_nidxslave);

        const auto iter_maxdist = s_elem.seek("maxdist");

        if (iter_maxdist == s_elem.end()) {
                throw std::runtime_error("elements.sfncon{4|6}.maxdist not defined");
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
                        throw std::runtime_error("elements.sfncon{4|6}.master must be an integer array");
                }

                const int32NDArray nidxslave = ov_nidxslave(l).int32_array_value();

                if (error_state) {
                        throw std::runtime_error("elements.sfncon{4|6}.slave must be an integer array");
                }

                ColumnVector maxdist = ov_maxdist(l).column_vector_value();

                if (error_state) {
                        throw std::runtime_error("elements.sfncon{4|6}.maxdist must be a column vector");
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
                default:
                        FEM_ASSERT(false);
                }

                if (nidxmaster.ndims() != 2 || nidxmaster.rows() < 1 || nidxmaster.columns() != iNumNodesElem) {
                        throw std::runtime_error("elements.sfncon{4|6}.master must be an nx{4|6} array");
                }

                if (nidxslave.ndims() != 2 || nidxslave.rows() < 1 || nidxslave.columns() != 1) {
                        throw std::runtime_error("elements.sfncon{4|6}.slave must be an nx1 array");
                }

                if (maxdist.rows() == 1 && nidxslave.rows() > 1) {
                        const double maxdistval = maxdist(0);
                        maxdist.resize(nidxslave.rows(), maxdistval);
                }

                if (maxdist.rows() != nidxslave.rows()) {
                        throw std::runtime_error("elements.sfncon{4|6}.maxdist must have "
                                                 "the same dimensions like elements.sfncon{4|6}.slave");
                }

                for (octave_idx_type i = 0; i < nidxmaster.rows(); ++i) {
                        for (octave_idx_type j = 0; j < nidxmaster.columns(); ++j) {
                                if (nidxmaster(i, j).value() < 1 || nidxmaster(i, j).value() > nodes.rows()) {
                                        throw std::runtime_error("elements.sfncon{4|6}.master: node index out of range");
                                }
                        }
                }

                for (octave_idx_type i = 0; i < nidxslave.rows(); ++i) {
                        if (nidxslave(i).value() < 1 || nidxslave(i).value() > nodes.rows()) {
                                throw std::runtime_error("elements.sfncon{4|6}.slave: node index out of range");
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
// PKG_ADD: autoload("FEM_SCA_STRESS_VMIS", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_CT_FIXED", "__mboct_fem_pkg__.oct");
// PKG_ADD: autoload("FEM_CT_SLIDING", "__mboct_fem_pkg__.oct");

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
// PKG_DEL: autoload("FEM_SCA_STRESS_VMIS", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_CT_FIXED", "__mboct_fem_pkg__.oct", "remove");
// PKG_DEL: autoload("FEM_CT_SLIDING", "__mboct_fem_pkg__.oct", "remove");


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

        if (locked_dof.ndims() != 2 ||
            locked_dof.columns() != nodes.columns() ||
            locked_dof.rows() != nodes.rows()) {
                error("size of load_case.locked_dof does not match size of mesh.nodes");
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
                case ElementTypes::ELEM_TET10: {
                        const int32NDArray elnodes = m_elements.contents(iter_elem_type).int32_array_value();

                        if (error_state) {
                                return retval;
                        }

                        for (octave_idx_type j = 0; j < elnodes.rows(); ++j) {
                                for (octave_idx_type k = 0; k < elnodes.columns(); ++k) {
                                        const octave_idx_type idxnode = elnodes(j, k).value();

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
                                                dof_in_use(idxnode - 1, l) = true;
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

        int32NDArray ndof(dim_vector(nodes.rows(), nodes.columns()), -1);

        for (octave_idx_type i = 0; i < ndof.rows(); ++i) {
                for (octave_idx_type j = 0; j < ndof.columns(); ++j) {
                        if (dof_in_use(i, j)) {
                                ndof(i, j) = locked_dof(i, j) ? 0 : ++icurrdof;
                        }
                }
        }

        octave_map dof_map;

        dof_map.assign("ndof", octave_value(ndof));

        int32NDArray idx_node(dim_vector(icurrdof, 1), -1);
        octave_idx_type icurrndof = 0;

        for (octave_idx_type i = 0; i < ndof.rows(); ++i) {
                for (octave_idx_type j = 0; j < ndof.columns(); ++j) {
                        if (ndof(i, j).value() > 0) {
                                idx_node(icurrndof++) = ndof(i, j);
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
                        case ElementTypes::ELEM_SFNCON6: {
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
                                case ElementTypes::ELEM_SFNCON6: {
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
                                                                edof[CS_JOINT].dof(edof[CS_JOINT].elem_idx, l) = ++icurrdof;
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
                                                        edof[CS_RBE3].dof(j, l) = ++icurrdof;
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
                                const octave_idx_type idxedof = edof[k].dof(i, j).value();

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

                for (octave_idx_type k = ElementTypes::ELEM_SFNCON4; k <= ElementTypes::ELEM_SFNCON6; ++k) {
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
                error(err.what());
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

                if (ndof.rows() != nodes.rows() || ndof.columns() != nodes.columns()) {
                        throw std::runtime_error("shape of dof_map.ndof does not match shape of nodes in argument dof_map");
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

                const auto iterC = material_data.seek("C");
                const auto iterRho = material_data.seek("rho");
                const auto iterAlpha = material_data.seek("alpha");
                const auto iterBeta = material_data.seek("beta");

                if (iterC == material_data.end()) {
                        throw std::runtime_error("field \"C\" not found in mesh.material_data in argument mesh");
                }

                if (iterRho == material_data.end()) {
                        throw std::runtime_error("field \"rho\" not found in mesh.material_data in argument mesh");
                }

                const Cell cellC = material_data.contents(iterC);
                const Cell cellRho = material_data.contents(iterRho);
                const Cell cellAlpha = iterAlpha != material_data.end() ? material_data.contents(iterAlpha) : Cell();
                const Cell cellBeta = iterBeta != material_data.end() ? material_data.contents(iterBeta) : Cell();

                vector<Material> rgMaterials;

                rgMaterials.reserve(material_data.numel());

                Matrix C(6, 6);

                for (octave_idx_type i = 0; i < material_data.numel(); ++i) {
                        C = cellC(i).matrix_value();

                        const double rho = cellRho(i).scalar_value();

                        if (error_state) {
                                throw std::runtime_error("mesh.material_data.rho is not a valid scalar in argument mesh");
                        }

                        if (C.rows() != 6 || C.columns() != 6) {
                                throw std::runtime_error("size of constitutive matrix mesh.material_data.C is not valid in argument mesh");
                        }

                        const double alpha = iterAlpha != material_data.end() ? cellAlpha(i).scalar_value() : 0.;

                        if (error_state) {
                                throw std::runtime_error("mesh.material_data.alpha is not a valid scalar in argument mesh");
                        }

                        const double beta = iterBeta != material_data.end() ? cellBeta(i).scalar_value() : 0.;

                        if (error_state) {
                                throw std::runtime_error("mesh.material_data.beta is not a valid scalar in argument mesh");
                        }

                        rgMaterials.emplace_back(C, rho, alpha, beta);
                }

                DofMap oDof(ndof, edof, inumdof);

                array<bool, ElementTypes::iGetNumTypes()> rgElemUse;

                std::fill(std::begin(rgElemUse), std::end(rgElemUse), false);

                for (octave_idx_type i = 0; i < matrix_type.numel(); ++i) {
                        switch (matrix_type(i).value()) {
                        case Element::MAT_STIFFNESS:
                        case Element::MAT_STIFFNESS_SYM:
                        case Element::MAT_STIFFNESS_SYM_L:
                                rgElemUse[ElementTypes::ELEM_RBE3] = true;
                                rgElemUse[ElementTypes::ELEM_JOINT] = true;
                                rgElemUse[ElementTypes::ELEM_SFNCON4] = true;
                                rgElemUse[ElementTypes::ELEM_SFNCON6] = true;
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
                        case Element::VEC_STRESS_CAUCH:
                        case Element::SCA_STRESS_VMIS:
                                rgElemUse[ElementTypes::ELEM_ISO8] = true;
                                rgElemUse[ElementTypes::ELEM_TET10] = true;
                                break;

                        case Element::VEC_LOAD_CONSISTENT:
                        case Element::VEC_LOAD_LUMPED:
                                if (load_case.numel() == 0) {
                                        throw std::runtime_error("missing argument load_case for matrix_type == FEM_VEC_LOAD_*");
                                }

                                rgElemUse[ElementTypes::ELEM_PRESSURE_ISO4] = true;
                                rgElemUse[ElementTypes::ELEM_PRESSURE_TRIA6] = true;
                                rgElemUse[ElementTypes::ELEM_STRUCT_FORCE] = true;
                                rgElemUse[ElementTypes::ELEM_JOINT] = true;
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
                        case ElementTypes::ELEM_TET10: {
                                const auto iter_elem = elements.seek(oElemType.name);

                                if (iter_elem == elements.end()) {
                                        continue;
                                }

                                int32NDArray elem_nodes = elements.contents(iter_elem).int32_array_value();

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
                                        rgElemBlocks.emplace_back(new ElementBlock<Iso8>(oElemType.type, elem_nodes, nodes, 3, elem_mat, rgMaterials));
                                        break;

                                case ElementTypes::ELEM_TET10:
                                        rgElemBlocks.emplace_back(new ElementBlock<Tet10>(oElemType.type, elem_nodes, nodes, 3, elem_mat, rgMaterials));
                                        break;

                                default:
                                        throw std::runtime_error("invalid element type");
                                }
                        } break;
                        case ElementTypes::ELEM_RBE3:
                        case ElementTypes::ELEM_JOINT: {
                                const auto iter_elem = elements.seek(oElemType.name);

                                if (iter_elem == elements.end()) {
                                        continue;
                                }

                                const octave_map s_elem(elements.contents(iter_elem).map_value());

                                if (error_state) {
                                        throw std::runtime_error("mesh.elements.rbe3 must be an struct array in argument mesh");
                                }

                                const auto iter_nodes = s_elem.seek("nodes");

                                if (iter_nodes == s_elem.end()) {
                                        throw std::runtime_error("missing field mesh.elements.rbe3.nodes in argument mesh");
                                }

                                const Cell ov_nodes(s_elem.contents(iter_nodes));

                                FEM_ASSERT(ov_nodes.numel() == s_elem.numel());

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
                                                throw std::runtime_error("mesh.elements.rbe3.nodes must be an array of integers in argument mesh");
                                        }

                                        if (elem_nodes.columns() < oElemType.min_nodes) {
                                                throw std::runtime_error(std::string("invalid number of nodes for element type ") + oElemType.name + " in argument mesh");
                                        }

                                        if (elem_nodes.rows() != 1) {
                                                throw std::runtime_error(std::string("invalid number of rows in node matrix for element type ") + oElemType.name + " in argument mesh");
                                        }

                                        for (octave_idx_type j = 0; j < elem_nodes.columns(); ++j) {
                                                octave_idx_type inode = elem_nodes(j);

                                                if (inode < 1 || inode > nodes.rows()) {
                                                        throw std::runtime_error(std::string("invalid node index for element type ") + oElemType.name + " in argument mesh");
                                                }
                                        }

                                        Matrix X(6, elem_nodes.columns());

                                        for (octave_idx_type j = 0; j < elem_nodes.columns(); ++j) {
                                                for (octave_idx_type k = 0; k < X.rows(); ++k) {
                                                        X(k, j) = nodes(elem_nodes(j).value() - 1, k);
                                                }
                                        }

                                        switch (oElemType.type) {
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
                                                                        U(l, k) = Uk(l);
                                                                }
                                                        }
                                                }

                                                pElem->Insert<ElemJoint>(++dofelemid[oElemType.dof_type], X, nullptr, elem_nodes, C, U);
                                        } break;
                                        default:
                                                FEM_ASSERT(false);
                                        }
                                }

                                if (dofelemid[oElemType.dof_type] > edof[oElemType.dof_type].rows()) {
                                        throw std::runtime_error("dof_map.edof is not consistent with elements in argument dof_map");
                                }

                                rgElemBlocks.emplace_back(std::move(pElem));
                        } break;
                        case ElementTypes::ELEM_SFNCON4:
                        case ElementTypes::ELEM_SFNCON6: {
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
                        case ElementTypes::ELEM_PRESSURE_TRIA6:
                                InsertPressureElem<PressureLoadTria6>(oElemType.type, nodes, load_case, oElemType.name, oElemType.max_nodes, rgElemBlocks);
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

                                                                                X(m, l) = nodes(inode, m);
                                                                        }
                                                                }

                                                                pElem->Insert(i + 1, X, nullptr, loaded_nodes, i + 1, loads);
                                                        }
                                                }

                                                rgElemBlocks.emplace_back(std::move(pElem));
                                        }
                                }
                        } break;
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
                        const Element::MatrixType eMatType = static_cast<Element::MatrixType>(matrix_type(j).value());
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
                        const Element::MatrixType eMatType = static_cast<Element::MatrixType>(matrix_type(i).value());

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
                        case Element::MAT_ACCEL_LOAD: {
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
                                                        mat(i, j) = mat(j, i);
                                                }
                                        }
                                        break;

                                default:
                                        break;
                                }

                                retval.append(mat);
                        } break;
                        case Element::VEC_STRESS_CAUCH:
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

                                octave_scalar_map s_tau, s_taum, s_vmis;

                                for (octave_idx_type j = 0; j < ElementTypes::iGetNumTypes(); ++j) {
                                        const ElementTypes::TypeInfo& oElemType = ElementTypes::GetType(j);

                                        if (!rgElemUse[oElemType.type]) {
                                                continue;
                                        }

                                        switch (oElemType.type) {
                                        case ElementTypes::ELEM_ISO8:
                                        case ElementTypes::ELEM_TET10: {
                                                const auto iter_elem = elements.seek(oElemType.name);

                                                if (iter_elem == elements.end()) {
                                                        continue;
                                                }

                                                const int32NDArray elem_nodes = elements.contents(iter_elem).int32_array_value();

                                                NDArray tau(dim_vector(elem_nodes.rows(), elem_nodes.columns(), iNumStress, iNumLoads), 0.);

                                                for (auto k = rgElemBlocks.cbegin(); k != rgElemBlocks.cend(); ++k) {
                                                        if ((*k)->GetElementType() == oElemType.type) {
                                                                (*k)->PostProcElem(tau, Element::VEC_STRESS_CAUCH, sol_U);
                                                        }
                                                }

                                                int32NDArray itaun(dim_vector(nodes.rows(), 1), 0);

                                                for (octave_idx_type k = 0; k < tau.dim1(); ++k) {
                                                        for (octave_idx_type l = 0; l < tau.dim2(); ++l) {
                                                                const octave_idx_type inode = elem_nodes(k, l).value() - 1;
                                                                itaun(inode) += 1;
                                                        }
                                                }

                                                NDArray taun(dim_vector(nodes.rows(), iNumStress, iNumLoads), 0.);

                                                for (octave_idx_type n = 0; n < iNumLoads; ++n) {
                                                        for (octave_idx_type m = 0; m < iNumStress; ++m) {
                                                                for (octave_idx_type k = 0; k < tau.dim1(); ++k) {
                                                                        for (octave_idx_type l = 0; l < tau.dim2(); ++l) {
                                                                                const octave_idx_type inode = elem_nodes(k, l).value() - 1;
                                                                                taun(inode, m, n) += tau(k, l, m + n * iNumStress) / itaun(inode).value();
                                                                        }
                                                                }
                                                        }
                                                }

                                                NDArray taum(tau.dims());

                                                for (octave_idx_type n = 0; n < iNumLoads; ++n) {
                                                        for (octave_idx_type m = 0; m < iNumStress; ++m) {
                                                                for (octave_idx_type k = 0; k < tau.dim1(); ++k) {
                                                                        for (octave_idx_type l = 0; l < tau.dim2(); ++l) {
                                                                                const octave_idx_type inode = elem_nodes(k, l).value() - 1;
                                                                                taum(k, l, m + n * iNumStress) = taun(inode, m, n);
                                                                        }
                                                                }
                                                        }
                                                }

                                                if (eMatType == Element::SCA_STRESS_VMIS) {
                                                        NDArray vmis(dim_vector(elem_nodes.rows(), elem_nodes.columns(), iNumLoads));

                                                        for (octave_idx_type n = 0; n < iNumLoads; ++n) {
                                                                for (octave_idx_type l = 0; l < tau.dim2(); ++l) {
                                                                        for (octave_idx_type k = 0; k < tau.dim1(); ++k) {
                                                                                const octave_idx_type ioffset = n * iNumStress;

                                                                                const double tauxx = taum(k, l, ioffset);
                                                                                const double tauyy = taum(k, l, ioffset + 1);
                                                                                const double tauzz = taum(k, l, ioffset + 2);
                                                                                const double tauxy = taum(k, l, ioffset + 3);
                                                                                const double tauyz = taum(k, l, ioffset + 4);
                                                                                const double tauzx = taum(k, l, ioffset + 5);

                                                                                vmis(k, l, n) = sqrt(tauxx * tauxx + tauyy * tauyy + tauzz * tauzz
                                                                                                     - (tauxx * tauyy + tauyy * tauzz + tauxx * tauzz)
                                                                                                     + 3. * (tauxy * tauxy + tauyz * tauyz + tauzx * tauzx));

                                                                        }
                                                                }
                                                        }

                                                        s_vmis.assign(oElemType.name, vmis);
                                                } else {
                                                        s_tau.assign(oElemType.name, tau);
                                                        s_taum.assign(oElemType.name, taum);
                                                }
                                        } break;

                                        default:
                                                break;
                                        }
                                }

                                octave_scalar_map stress;

                                if (eMatType == Element::SCA_STRESS_VMIS) {
                                        stress.assign("vmis", s_vmis);
                                } else {
                                        stress.assign("tau", s_tau);
                                        stress.assign("taum", s_taum);
                                }

                                retval.append(stress);
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
                        warning_with_id("mboct-fem-pkg:invalid-mesh", "Jacobian is singular or negative det(J)=%g", detJmin);
                }
        } catch (const std::exception& err) {
                error(err.what());
                return retval;
        }

        return retval;
}

#define DEFINE_GLOBAL_CONSTANT(NAMESPACE, CONST, DESCRIPTION)           \
        DEFUN_DLD(FEM_##CONST, args, nargout,                           \
                  "-*- texinfo -*-\n"                                   \
                  "@deftypefn {} @var{id} = FEM_" #CONST  "()\n"        \
                  DESCRIPTION "\n"                                      \
                  "@end deftypefn\n")                                   \
                                                                        \
        {                                                               \
                return octave_value(NAMESPACE::CONST);                  \
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
DEFINE_GLOBAL_CONSTANT(Element, VEC_STRESS_CAUCH, "Cauchy stress")
DEFINE_GLOBAL_CONSTANT(SurfToNodeConstrBase, CT_FIXED, "build constraints in all three directions in space")
DEFINE_GLOBAL_CONSTANT(SurfToNodeConstrBase, CT_SLIDING, "build only one constraint normal to the surface")