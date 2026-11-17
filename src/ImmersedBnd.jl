abstract type AbstractBoundary end
abstract type AbstractBoundaryFunction end

#This is reserved for a future version of the code. 
abstract type AbstractIntegralOperator end

#Lets define an Eulerian Grid
mutable struct PeriodicLagBnd <: AbstractBoundary
	#Size of the immersed boundary
	N::Int
	#Location of immersed boundary points, as well as their forces
	bnd_loc::Vector{Point2}
	#Reference location. Possibly useful, unclear. 
	ref_bnd_loc::Vector{Point2}
	#Current force density. 
	stretch_force::Vector{Point2}
	# connect_force::Matrix{Float64} #This one won't be used until we have an immersed mesh. 
	#Arrays for constants that define constitutive laws for the above forces
	stretch_kp::Vector{Float64}
	stretch_km::Vector{Float64}
	gammap::Vector{Float64}
	gammam::Vector{Float64}
	
	#Arrrays which result from finite differences to approximate tangent vectors
	#REFERENCE CONFIGURATION
	ref_taup::Vector{Point2}
	ref_taum::Vector{Point2}
	ref_tau0::Vector{Point2}
	ref_norm0::Vector{Point2}
	#Arrays from finite differences to approximate arc length
	ref_dsp::Vector{Float64}
	ref_dsm::Vector{Float64}
	ref_ds0::Vector{Float64}

	#CURRENT CONFIGURATION
	#Arrrays which result from finite differences to approximate tangent vectors
	taup::Vector{Point2}
	taum::Vector{Point2}
	tau0::Vector{Point2}
	norm0::Vector{Point2}
	#Arrays from finite differences to approximate arc length
	dsp::Vector{Float64}
	dsm::Vector{Float64}
	ds0::Vector{Float64}

	#Arrays which are used to flag and define points of the immersed network to which the bnd is connected
	# flag::Vector{Int}  #This one won't be used until we have an immersed mesh. 
	# connect_buddy{Int} #This one won't be used until we have an immersed mesh. 

end


#We can make a periodic lagrangian boundary as long as we get two vectors of location points. 
function PeriodicLagBnd(X::Vector{<:Real},Y::Vector{<:Real})
	if ~(size(X) == size(Y))
		throw(ArgumentError("Size of inputs does not match"));
	end
	N = length(X);
	location = Point2.(X,Y);
	
	#Need an array to hold the forces. Initialize at zero. 
	stretch_force = zeros(Point2, N);

	kp = ones(N);
	km = ones(N);
	gammap = ones(N);
	gammam = ones(N);

	#Tangents, arc lengths and (inward) normals come from the helper in LagGeometry.jl. 
	#THE NORMAL CALCULATION ASSUMES THAT WE HAVE A COUNTER-CLOCK-WISE
	#ORIENTATION OF A CLOSED CURVE. 
	taup,taum,tau0,norm0,dsp,dsm,ds0 = PeriodicCurveGeometry(location);

	#Set reference geometry equal to current. 
	ref_bnd_loc = copy(location);
	ref_taup = copy(taup);
	ref_taum = copy(taum);
	ref_tau0 = copy(tau0);
	ref_norm0 = copy(norm0);
	ref_dsp = copy(dsp);
	ref_dsm = copy(dsm);
	ref_ds0 = copy(ds0);

	bnd::PeriodicLagBnd = PeriodicLagBnd(N,location,ref_bnd_loc,stretch_force,kp,km,gammap,gammam,
										ref_taup,ref_taum,ref_tau0,ref_norm0,ref_dsp,ref_dsm,
										ref_ds0,taup,taum,tau0,norm0,dsp,dsm,ds0);
	return bnd

end

#This is a struct which holds data defined on the boundary. First scalar
mutable struct ScalarBndData <: AbstractBoundaryFunction
	U::Vector{Float64}
	bnd::AbstractBoundary

	function ScalarBndData(u::Vector{Float64},mybnd::AbstractBoundary)
		if ~(length(u) == mybnd.N)
			throw(ArgumentError("Size of input does not match boundary size"));
		else
			return new(u,mybnd)
		end
	end
end

#And now vector
mutable struct VectorBndData <: AbstractBoundaryFunction
	U::Vector{Float64}
	V::Vector{Float64}
	bnd::AbstractBoundary

	function VectorBndData(u::Vector{Float64},v::Vector{Float64},mybnd::AbstractBoundary)
		if ~(length(u) == mybnd.N)
			throw(ArgumentError("Size of first input does not match boundary size"));
		elseif ~(length(v) == mybnd.N)
			throw(ArgumentError("Size of second input does not match boundary size"));
		else
			return new(u,v,mybnd)
		end
	end
end




############################################
#  Now for some spread and inter operators #
############################################

function ScalarBndSpread(data::Vector{Float64},mybnd::AbstractBoundary,mygrid::AbstractGrid)
	if ~(length(data) == mybnd.N)
		throw(ArgumentError("Size of data does not match boundry"))
	end

	M = mybnd.N;
	Nx = mygrid.Nx;
	Ny = mygrid.Ny;

	output = zeros(Nx,Ny);

	for i = 1:M

		ix, xweights = DeltaWeights(mybnd.bnd_loc[i][1], mygrid.Lx, mygrid.dx, Nx)
		iy, yweights = DeltaWeights(mybnd.bnd_loc[i][2], mygrid.Ly, mygrid.dy, Ny)

		for (j,indx) in pairs(ix)
			for (k,indy) in pairs(iy)
				output[indx,indy] = output[indx,indy] + data[i]*xweights[j]*yweights[k]*mybnd.ref_ds0[i]
			end
		end

	end
	return output
end


function ScalarBndInterp(data::Matrix{Float64},mybnd::AbstractBoundary,mygrid::AbstractGrid)
	if ~(size(data) == (mygrid.Nx,mygrid.Ny))
		throw(ArgumentError("Size of data does not match grid"))
	end

	M = mybnd.N;
	Nx = mygrid.Nx;
	Ny = mygrid.Ny;

	output = zeros(M);

	for i = 1:M

		ix, xweights = DeltaWeights(mybnd.bnd_loc[i][1], mygrid.Lx, mygrid.dx, Nx)
		iy, yweights = DeltaWeights(mybnd.bnd_loc[i][2], mygrid.Ly, mygrid.dy, Ny)

		for (j,indx) in pairs(ix)
			for (k,indy) in pairs(iy)
				output[i] = output[i] + data[indx,indy]*xweights[j]*yweights[k]*mygrid.dx*mygrid.dy
			end
		end

	end
	return output
end




function PeskinDelta(X::AbstractVector{<:Real})
	if length(X) != 4
		throw(ArgumentError("PeskinDelta expects the 4 stencil distances"))
	end
	delta = zeros(4)
	delta[1] = (5.0 + 2.0*X[1] - sqrt(-7. - 12.0*X[1] - 4.0*X[1]^2) )/8.0
	delta[2] = (3.0 + 2.0*X[2] + sqrt(1. - 4.0*X[2] - 4.0*X[2]^2) )/8.0
	delta[3] = (3.0 - 2.0*X[3] + sqrt(1. + 4.0*X[3] - 4.0*X[3]^2) )/8.0
	delta[4] = (5.0 - 2.0*X[4] - sqrt(-7. + 12.0*X[4] - 4.0*X[4]^2) )/8.0
	return delta
end



# Wrapped grid indices and delta weights for a point at position `pos` along
# one axis of a cell-centered periodic grid that starts at `lo` with spacing `h`
# and `N` cells. The four stencil points are ordered left to right.
function DeltaWeights(pos::Real, lo::Real, h::Real, N::Integer)
	if !(h > 0)
		throw(ArgumentError("grid spacing must be positive"))
	end
	if !(N>0)
		throw(ArgumentError("number of grid cells must be positive"))
	end
	if !(isfinite(pos))
		throw(ArgumentError("boundary point position is not finite: $pos"))
	end

	left = floor(Int, (pos - lo)/h + 0.5)
	idx = left-1:left+2
	centers = lo .+ h*(idx .- 0.5)
	weights = PeskinDelta((centers .- pos)/h)./h
	return mod1.(idx, N), weights
end



###########################################################
# Simple integration techniques will be necessary as well #
###########################################################
#On periodic grids simple quadrature rules should have an additional
#Order of accuracy. Should be good enough for me. 

function  BndIntegral(data::Vector{Float64},mybnd::AbstractBoundary; reference::Bool=true)
	if ~(length(data) == mybnd.N)
		throw(ArgumentError("Size of data does not match boundry"))
	end

	if reference
		ds = mybnd.ref_ds0;
	else
		ds = mybnd.ds0;
	end

	summand = @. data*ds;
	output = sum(summand);
	return output
end