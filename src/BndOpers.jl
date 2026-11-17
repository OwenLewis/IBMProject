# A collection of functions associated with immersed boundary data types. 
# This is where force calculations, spread and interp, as well as 
import HDF5


# This is a simple function to read a network in from an HDF5 file
function ReadBnd(filename::String)
	fid = HDF5.h5open(filename)

	X = vec(HDF5.read(fid,"X"));
	Y = vec(HDF5.read(fid,"Y"));

	mybnd = PeriodicLagBnd(X,Y);
	close(fid)

	return mybnd
end


# These functions update all of the geometry characteristics of a periodic boundary. 
# The assumption is that all of the point coords have been updated so that 
# They are no longer equal to the ref_bnd_loc
function DefUpdate!(mybnd::PeriodicLagBnd)
	#All of the geometry comes from the helper in LagGeometry.jl
	mybnd.taup, mybnd.taum, mybnd.tau0, mybnd.norm0, mybnd.dsp, mybnd.dsm, mybnd.ds0 = PeriodicCurveGeometry(mybnd.bnd_loc);
	return nothing
end


# This function moves a boundary by a given displacement field. 
# It comes in three types

#Two distinct arrays with x&y displacements
function MoveBnd!(mybnd::PeriodicLagBnd,x::Vector{<:Real},y::Vector{<:Real})
	# First lets have some sanity checks to make sure that 
	if ~(size(x) == size(y))
		throw(ArgumentError("Size of displacements do not match"));
	end
	if ~(length(x) == mybnd.N)
		throw(ArgumentError("Size of displacements does not match bnd points"));
	end

	# Update the coordinates of all the veritices 
	mybnd.bnd_loc .+= Point2.(x,y);

	# Then update the mesh
	DefUpdate!(mybnd);
end

#An array of displacement vectors. 
function MoveBnd!(mybnd::PeriodicLagBnd,X::Vector{<: Vector{<:Real}})
	# First lets have some sanity checks to make sure that 
	if ~(length(X[1]) == 2)
		throw(ArgumentError("Size of displacements is not appropriate"));
	end
	if ~(length(X) == mybnd.N)
		throw(ArgumentError("Size of displacements does not match boundary points"));
	end

	# Update the coordinates of all the veritices 
	mybnd.bnd_loc .+= Point2.(X);
	# Then update the mesh
	DefUpdate!(mybnd);
end


#And the same for displacements that are already stored as Point2 (static 2-vectors)
function MoveBnd!(mybnd::PeriodicLagBnd,X::Vector{<:StaticArrays.SVector{2,<:Real}})
	if ~(length(X) == mybnd.N)
		throw(ArgumentError("Size of displacements does not match boundary points"));
	end

	mybnd.bnd_loc .+= Point2.(X);
	DefUpdate!(mybnd);
end


function ElasticForceCalc(mybnd::PeriodicLagBnd)
	#This calculates the elastic force in the PeriodicLagBnd due to its current configuration
	ForceDensity = 0;
	#JUST A PLACE-HOLDER FOR NOW!!
	return ForceDensity
end

