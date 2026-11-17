# Geometry on plain 2D points. Knows nothing about mesh/bnd types,
# so both the constructors (ImmersedNetwork.jl/ImmersedBnd.jl) and DefUpdate! 
# (NetworkOpers.jl/BndOpers.jl) can use it.
# Defines new type Point2 (vector in the plane) which other types can use. 


import LinearAlgebra
import StaticArrays

# One point in the plane. All point-valued data on the immersed boundary is stored as Vector{Point2}.
const Point2 = StaticArrays.SVector{2,Float64}

# Area of the triangle p1,p2,p3. The abs() is deliberate: the area is unsigned,
# so an inverted triangle is detected separately in DefUpdate! through J = det(Xmat)/det(Smat).
TriangleArea(p1,p2,p3) = abs(p1[1]*(p2[2] - p3[2]) +
                             p2[1]*(p3[2] - p1[2]) +
                             p3[1]*(p1[2] - p2[2]))/2

TriangleCentroid(p1,p2,p3) = (p1 + p2 + p3)/3

# 2x2 matrix whose columns are the edge vectors p2-p1 and p3-p1.
# Gives Xmat when called with current coordinates and Smat with reference coordinates.
TriangleEdgeMatrix(p1,p2,p3) = hcat(p2 .- p1, p3 .- p1)


# Simple Midpoint calculation of two points in the plane. Can also be used to average vectors. 
Midpoint(p1,p2) = (p1 .+ p2)./2

#Simple rotation (counterclockwise) of a vector in the plane by 90 degress
#Rotating the tangent gives the INWARD normal for a counter-clockwise curve
Rotate90(t) = Point2(-t[2], t[1])


# Tangent, arc-length and normal data for a closed periodic curve whose points are given as a Vector{Point2}.
function PeriodicCurveGeometry(loc::AbstractVector{<:StaticArrays.SVector{2,<:Real}})
    taup = circshift(loc,-1) .- loc          #forward difference, wraps around (periodicity)
    taum = loc .- circshift(loc,1)           #backward difference, wraps around
    tau0 = Midpoint.(taup,taum)

    dsp = LinearAlgebra.norm.(taup)
    dsm = LinearAlgebra.norm.(taum)
    ds0 = Midpoint.(dsp,dsm)

    #Rotate the tangent to get normal
    norm0 = Rotate90.(tau0)

    return taup, taum, tau0, norm0, dsp, dsm, ds0
end
