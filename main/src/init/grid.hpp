/*
 * MIT License
 *
 * Copyright (c) 2021 CSCS, ETH Zurich
 *               2021 University of Basel
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

/*! @file
 * @brief Initial placement of particles on a regular grid
 *
 * @author Sebastian Keller <sebastian.f.keller@gmail.com>
 */

#pragma once

#include "cstone/sfc/box.hpp"
#include "cstone/sfc/sfc.hpp"
#include "cstone/util/array.hpp"
#include "cstone/util/gsl-lite.hpp"

namespace sphexa
{

/*! @brief partition a range R into N segments
 *
 * @param R   the length of the range to be divided
 * @param i   the index of the segment to be computed
 * @param N   the total number of segments
 * @return    the start and end index (subrange of R) of the i-th segment
 */
auto partitionRange(size_t R, size_t i, size_t N)
{
    size_t s = R / N;
    size_t r = R % N;
    if (i < r)
    {
        size_t start = (s + 1) * i;
        size_t end   = start + s + 1;
        return std::make_tuple(start, end);
    }
    else
    {
        size_t start = (s + 1) * r + s * (i - r);
        size_t end   = start + s;
        return std::make_tuple(start, end);
    }
}

/*! @brief create regular cubic grid centered on (0,0,0), spanning [-r, r)^3, in the x,y,z arrays
 *
 * @tparam     Vector
 * @param[in]  r      half the cube side-length
 * @param[in]  side   number of particles along each dimension
 * @param[in]  first  index in [0, side^3] of first particle add to x,y,z
 * @param[in]  last   index in [first, side^3] of last particles to add to x,y,z
 * @param[out] x      output coordinates, length = last - first
 * @param[out] y
 * @param[out] z
 */
template<class Vector>
void regularGrid(double r, size_t side, size_t first, size_t last, Vector& x, Vector& y, Vector& z)
{
    double step = (2. * r) / side;

#pragma omp parallel for
    for (size_t i = first / (side * side); i < last / (side * side) + 1; ++i)
    {
        double lz = -r + (i * step);

        for (size_t j = 0; j < side; ++j)
        {
            double ly = -r + (j * step);

            for (size_t k = 0; k < side; ++k)
            {
                size_t lindex = (i * side * side) + (j * side) + k;

                if (first <= lindex && lindex < last)
                {
                    double lx = -r + (k * step);

                    z[lindex - first] = lz;
                    y[lindex - first] = ly;
                    x[lindex - first] = lx;
                }
            }
        }
    }
}


/*! @brief create regular cubic grid centered on (0,0,0) with an internal cube and an external cube with different rho
 *
 * @tparam     Vector
 * @param[in]  r         half the internal cube side-length
 * @param[in]  rDelta    half the external delta cube side-length
 * @param[in]  stepInt   step in the internal part
 * @param[in]  stepExt   step in the external part
 * @param[in]  cubeSide  number of particles in the internal cube along each dimension
 * @param[in]  deltaSide number of particles in the internal cube along each dimension
 * @param[in]  first     index in [0, side^3] of first particle add to x,y,z
 * @param[in]  last      index in [first, side^3] of last particles to add to x,y,z
 * @param[out] x         output coordinates, length = last - first
 * @param[out] y
 * @param[out] z
 */
template<class Vector>
void internalCubeGrid(double r, double rDelta, double stepInt, double stepExt, size_t cubeSide, size_t deltaSide,
        size_t first, size_t last, Vector& x, Vector& y, Vector& z)
{
    size_t side  = cubeSide + (2. * deltaSide);

#pragma omp parallel for
    for (size_t i = first / (side * side); i < last / (side * side) + 1; ++i)
    {
        double lz;
        if (i < deltaSide)
        {
            lz = -(r + rDelta) + (i * stepExt);                          // i < -r
        }
        else if (i >= deltaSide && i <= deltaSide)
        {
            lz = -r + ((i-deltaSide) * stepInt);                         // r(i) >= -r && r(i) <= +r
        }
        else
        {
            lz = r + ((i-deltaSide-cubeSide) * stepExt);                 // r(i) > +r
        }

        for (size_t j = 0; j < side; ++j)
        {
            double ly;
            if (j < deltaSide)
            {
                ly = -(r + rDelta) + (j * stepExt);                      // j < -r
            }
            else if (j >= deltaSide && j <= deltaSide)
            {
                ly = -r + ((j-deltaSide) * stepInt);                     // r(j) >= -r && r(j) <= +r
            }
            else
            {
                ly = r + ((j-deltaSide-cubeSide) * stepExt);             // r(j) > +r
            }

            for (size_t k = 0; k < side; ++k)
            {
                size_t lindex = (i * side * side) + (j * side) + k;

                if (first <= lindex && lindex < last)
                {
                    double lx;
                    if (k < deltaSide)
                    {
                        lx = -(r + rDelta) + (k * stepExt);              // k < -r
                    }
                    else if (k >= deltaSide && k <= deltaSide)
                    {
                        lx = -r + ((k-deltaSide) * stepInt);             // r(k) >= -r && r(k) <= +r
                    }
                    else
                    {
                        lx = r + ((k-deltaSide-cubeSide) * stepExt);     // r(k) > +r
                    }


                    z[lindex - first] = lz;
                    y[lindex - first] = ly;
                    x[lindex - first] = lx;
                }
            }
        }
    }
}

/*! @brief intersection of a box with a regular grid
 *
 * @tparam T   float or double
 * @param  a   sub box of the unit-cube
 * @param  m   number of grid-segments per dimension
 * @return     two integer triples marking which grid cells intersected with @p a
 */
template<class T>
auto gridIntersection(const cstone::FBox<T>& a, int m)
{
    auto                l = util::array<T, 3>{a.xmin(), a.ymin(), a.zmin()} * T(m);
    util::array<int, 3> lowerIdx{int(l[0]), int(l[1]), int(l[2])};

    auto                u = util::array<T, 3>{a.xmax(), a.ymax(), a.zmax()} * T(m);
    util::array<int, 3> upperIdx{int(std::ceil(u[0])), int(std::ceil(u[1])), int(std::ceil(u[2]))};

    return std::make_tuple(lowerIdx, upperIdx);
}

/*! @brief scale coordinates from a [0,1] unit block into a global coordinate box
 *
 * @tparam T             float or double
 * @param  uX            3D coordinate in [0,1]^3
 * @param  gridIdx       3D integer coordinate in [0,m-1]^3
 * @param  m             multiplicity of the block-grid
 * @param  globalBox     the global coordinate bounds
 * @return               mapped coordinate into the specified grid index
 *
 * Here we have the globalBox decomposed into an m x m x m grid. Each of those m^3 components
 * can be indexed with @p gridIdx. We assume that @p uX comes from a template (glass) block with
 * coordinates in [0,1]. With this function we can map a template block into one of the m^3 grid points
 * of the global box.
 */
template<class T>
cstone::Vec3<T> scaleBlockToGlobal(cstone::Vec3<T> uX, cstone::Vec3<int> gridIdx, int m,
                                   const cstone::Box<T>& globalBox)
{

    cstone::Vec3<T> blockOrigin{gridIdx[0] * globalBox.lx(), gridIdx[1] * globalBox.ly(), gridIdx[2] * globalBox.lz()};
    blockOrigin /= T(m);

    cstone::Vec3<T> globalOrigin{globalBox.xmin(), globalBox.ymin(), globalBox.zmin()};

    auto gX = uX;
    gX[0] *= globalBox.lx();
    gX[1] *= globalBox.ly();
    gX[2] *= globalBox.lz();
    gX /= m;

    gX += globalOrigin + blockOrigin;

    return gX;
}

/*! @brief extract (push into vector) coordinates of the virtual global box that are contained in @p selectBox
 *
 * @tparam T
 * @tparam Vector
 * @param[in]  selectBox    a sub box of @p globalBox
 * @param[in]  globalBox    global coordinate bounding box
 * @param[in]  gridIdx      3D integer coordinate in [0,m-1]^3
 * @param[in]  m            multiplicity of the global grid
 * @param[in]  xBlock       x-coords of the template block in [0,1]
 * @param[in]  yBlock       y-coords of the template block in [0,1]
 * @param[in]  zBlock       z-coords of the template block in [0,1]
 * @param[out] x            output x-coords that lie in @p selectBox
 * @param[out] y            output y-coords that lie in @p selectBox
 * @param[out] z            output z-coords that lie in @p selectBox
 */
template<class T, class Vector>
void extractBlock(const cstone::FBox<T>& selectBox, const cstone::Box<T>& globalBox, cstone::Vec3<int> gridIdx, int m,
                  gsl::span<const T> xBlock, gsl::span<const T> yBlock, gsl::span<const T> zBlock, Vector& x, Vector& y,
                  Vector& z)
{
    for (size_t i = 0; i < xBlock.size(); ++i)
    {
        auto sX = scaleBlockToGlobal({xBlock[i], yBlock[i], zBlock[i]}, gridIdx, m, globalBox);

        bool select = (selectBox.xmin() <= sX[0] && sX[0] < selectBox.xmax()) &&
                      (selectBox.ymin() <= sX[1] && sX[1] < selectBox.ymax()) &&
                      (selectBox.zmin() <= sX[2] && sX[2] < selectBox.zmax());
        if (select)
        {
            x.push_back(sX[0]);
            y.push_back(sX[1]);
            z.push_back(sX[2]);
        }
    }
}

/*! @brief Construct a part of a coordinate box obtained by stacking a template block m-times in each dimension
 *
 * @tparam T
 * @tparam KeyType
 * @tparam Vector
 * @param[in]  keyStart     SFC key start
 * @param[in]  keyEnd       SFC key end
 * @param[in]  globalBox    global coordinate bounding box
 * @param[in]  multiplicity multiplicity of the global grid
 * @param[in]  xBlock       x-coords of the template block in [0,1]
 * @param[in]  yBlock       y-coords of the template block in [0,1]
 * @param[in]  zBlock       z-coords of the template block in [0,1]
 * @param[out] x            output x-coords with SFC keys in @p [keyStart:keyEnd]
 * @param[out] y            output y-coords with SFC keys in @p [keyStart:keyEnd]
 * @param[out] z            output z-coords with SFC keys in @p [keyStart:keyEnd]
 */
template<class T, class KeyType, class Vector>
void assembleCube(KeyType keyStart, KeyType keyEnd, const cstone::Box<T>& globalBox, int multiplicity,
                  gsl::span<const T> xBlock, gsl::span<const T> yBlock, gsl::span<const T> zBlock, Vector& x, Vector& y,
                  Vector& z)
{
    // span the assigned SFC range with valid octree cells
    int                  numCells = cstone::spanSfcRange(keyStart, keyEnd);
    std::vector<KeyType> cells(numCells + 1);
    cstone::spanSfcRange(keyStart, keyEnd, cells.data());
    cells.back() = keyEnd;

    // extract the volume of each cell from the virtual global glass block grid
    for (size_t i = 0; i < cstone::nNodes(cells); ++i)
    {
        auto iBox      = cstone::sfcIBox(cstone::sfcKey(cells[i]), cstone::sfcKey(cells[i + 1]));
        auto selectBox = cstone::createFpBox<KeyType>(iBox, globalBox);

        // determine which building blocks in the glass block grid the current selectBox intersects with
        auto [lowerIdx, upperIdx] =
            gridIntersection(cstone::createFpBox<KeyType>(iBox, cstone::Box<T>(0, 1)), multiplicity);
        for (int ix = lowerIdx[0]; ix < upperIdx[0]; ++ix)
            for (int iy = lowerIdx[1]; iy < upperIdx[1]; ++iy)
                for (int iz = lowerIdx[2]; iz < upperIdx[2]; ++iz)
                {
                    extractBlock<T>(selectBox, globalBox, {ix, iy, iz}, multiplicity, xBlock, yBlock, zBlock, x, y, z);
                }
    }
}

//! @brief Discard any particles outside a sphere with radius @p r centered on the origin
template<class T, class Vector>
void cutSphere(T radius, Vector& x, Vector& y, Vector& z)
{
    std::vector<int> particleSelection(x.size());

#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < x.size(); ++i)
    {
        T rp                 = std::sqrt(x[i] * x[i] + y[i] * y[i] + z[i] * z[i]);
        particleSelection[i] = (rp <= radius);
    }

    size_t numSelect = std::count(particleSelection.begin(), particleSelection.end(), 1);

    Vector xSphere, ySphere, zSphere;
    xSphere.reserve(numSelect);
    ySphere.reserve(numSelect);
    zSphere.reserve(numSelect);

    for (size_t i = 0; i < x.size(); ++i)
    {
        if (particleSelection[i])
        {
            xSphere.push_back(x[i]);
            ySphere.push_back(y[i]);
            zSphere.push_back(z[i]);
        }
    }

    swap(x, xSphere);
    swap(y, ySphere);
    swap(z, zSphere);
}

} // namespace sphexa
