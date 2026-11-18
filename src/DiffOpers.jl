# Code to be loaded into IBMProject
import FFTW
import LinearAlgebra

abstract type AbstractDifferentialOperator end


#Now we can define differential Operators which can act on grid data
#This one is simple and does not change the rank of the data
struct SimplePeriodicDifferentialOperator <: AbstractDifferentialOperator
	grid::PeriodicEulGrid
	Eigenvalues::Matrix{ComplexF64}
end

#This one is gradient-like. Applying it will act on scalars and produce vectors
#Inverting it will act on vectors and produce scalars
struct GradientPeriodicDifferentialOperator <: AbstractDifferentialOperator
	grid::PeriodicEulGrid
	EigenvaluesU::Matrix{ComplexF64}
	EigenvaluesV::Matrix{ComplexF64}
end

#This one is divergence-like. Applying it will act on vectors and produce scalars
#Inverting it will act on scalars and produce vectors
struct DivergencePeriodicDifferentialOperator <: AbstractDifferentialOperator
	grid::PeriodicEulGrid
	EigenvaluesU::Matrix{ComplexF64}
	EigenvaluesV::Matrix{ComplexF64}
end

struct InverseSimplePeriodicDifferentialOperator <: AbstractDifferentialOperator
	grid::PeriodicEulGrid
	InvEigenvalues::Matrix{ComplexF64}
	dropped::BitMatrix 
end


###############################################
#         Writing my operator algebra         #
###############################################
#These methods let me construct new operators from exiting ones via addition, 
#Scalar multiplication, and composition. 

## Simple Operators first
## Scalar Multiplication (from either side)
function Base.:*(a::Real,L::SimplePeriodicDifferentialOperator)
	return SimplePeriodicDifferentialOperator(L.grid,a.*L.Eigenvalues)
end
function Base.:*(L::SimplePeriodicDifferentialOperator,a::Real)
	return SimplePeriodicDifferentialOperator(L.grid,a.*L.Eigenvalues)
end

## Adddition
function Base.:+(A::SimplePeriodicDifferentialOperator,B::SimplePeriodicDifferentialOperator)
	CheckSameGrid(A,B,"added")
	return SimplePeriodicDifferentialOperator(A.grid,A.Eigenvalues.+B.Eigenvalues);
end

## Subtraction
function Base.:-(A::SimplePeriodicDifferentialOperator,B::SimplePeriodicDifferentialOperator)
	CheckSameGrid(A,B,"subtracted")
	return SimplePeriodicDifferentialOperator(A.grid,A.Eigenvalues.-B.Eigenvalues);
end
## Unary Minus
function Base.:-(A::SimplePeriodicDifferentialOperator) 
	return SimplePeriodicDifferentialOperator(A.grid, -A.Eigenvalues);
end

# Identity shifts, which could be useful for IBDL & IBSL. 
# Works for any REAL multiple of the identity
# This just avoids having to create an operator who's e-values are all a constant. 
function Base.:+(A::SimplePeriodicDifferentialOperator,J::LinearAlgebra.UniformScaling{<:Real})
	# λ is the field that stores the scalar multiple of the identity operator. 
	return SimplePeriodicDifferentialOperator(A.grid,A.Eigenvalues .+ J.λ)
end
function Base.:+(J::LinearAlgebra.UniformScaling{<:Real},A::SimplePeriodicDifferentialOperator)
	return A + J
end
function Base.:-(A::SimplePeriodicDifferentialOperator,J::LinearAlgebra.UniformScaling{<:Real})
	return A + (-J)
end
function Base.:-(J::LinearAlgebra.UniformScaling{<:Real},A::SimplePeriodicDifferentialOperator)
	return J + (-A)
end

# Now for operator composition. 
# There are only a few valid combinations

# Simple operators just get their e-values multiplied. 
function Base.:*(A::SimplePeriodicDifferentialOperator,B::SimplePeriodicDifferentialOperator)
	CheckSameGrid(A,B,"composed");
	return SimplePeriodicDifferentialOperator(A.grid,A.Eigenvalues.*B.Eigenvalues)
end

# This is a simple operator, composed with a gradient operator
# Result is a gradient operator
function Base.:*(G::GradientPeriodicDifferentialOperator,A::SimplePeriodicDifferentialOperator)
	CheckSameGrid(G,A,"composed")
	return GradientPeriodicDifferentialOperator(G.grid,G.EigenvaluesU.*A.Eigenvalues,G.EigenvaluesV.*A.Eigenvalues);
end

# Here we have a divergence operator, composed with a simple one
# Result is a divergence operator
function Base.:*(A::SimplePeriodicDifferentialOperator,D::DivergencePeriodicDifferentialOperator)
	CheckSameGrid(A,D,"composed");
	return DivergencePeriodicDifferentialOperator(A.grid,A.Eigenvalues.*D.EigenvaluesU,A.Eigenvalues.*D.EigenvaluesV);
end

# Here we have gradient, composed with divergence (grid laplacian)
# The result is a simple. 
function Base.:*(D::DivergencePeriodicDifferentialOperator,G::GradientPeriodicDifferentialOperator)
	CheckSameGrid(D,G,"composed")
	return SimplePeriodicDifferentialOperator(D.grid,D.EigenvaluesU.*G.EigenvaluesU+D.EigenvaluesV.*G.EigenvaluesV);
end

#Finally we have exponentiation 
function Base.:^(A::SimplePeriodicDifferentialOperator,n::Integer)
	if n < 1
		throw(ArgumentError("Power must be positive integer (use InvertOperator for inverse powers)."))
	end
	return SimplePeriodicDifferentialOperator(A.grid,A.Eigenvalues.^n);
end

###############################################
#     Now for some construction functions     #
###############################################


#This is to make a periodic laplacian of type PeriodicDifferentialOperator
function MakePeriodicLaplacian(mygrid::PeriodicEulGrid)
	ωx = FFTW.fftfreq(mygrid.Nx,mygrid.Nx/mygrid.L)
	ωy = FFTW.fftfreq(mygrid.Ny,mygrid.Ny/mygrid.H)
	compfreqx = im*ωx*(2*pi); #The eigenvalues of a single derivative in x&y
    compfreqy = im*ωy*(2*pi);
	# compfreqx = FirstDerivativeEigenvalues(mygrid.Nx,mygrid.L);
	# compfreqy = FirstDerivativeEigenvalues(mygrid.Ny,mygrid.H);

    ΩX = (compfreqx[i] for i=1:mygrid.Nx,j=1:mygrid.Ny); #put them into arrays the same size as the grid
    ΩY = (compfreqy[j] for i=1:mygrid.Nx,j=1:mygrid.Ny);

    applyeigs = (ΩX.^2 + ΩY.^2); #Eigenvalues of the laplacian
	Lap = SimplePeriodicDifferentialOperator(mygrid,applyeigs);
	return Lap
end


# Operator a*I - b*Laplacian. Its eigenvalues are a + b*|k|^2, since the Laplacian's are -|k|^2.
# Uses the true spectral Laplacian (which differs from D*G at Nyquist frequency on even grids).
function MakePeriodicHelmholtz(mygrid::PeriodicEulGrid, a::Real, b::Real)
	if a < 0 || b < 0
		throw(ArgumentError("a and b must be non-negative (otherwise the operator can be singular or indefinite)"));
	end
	if a == 0 && b == 0
		throw(ArgumentError("a and b cannot both be zero"));
	end
	Lap = MakePeriodicLaplacian(mygrid);
	return a*LinearAlgebra.I - b*Lap;
end

# Now the operators that change the rank of data. They don't technically store different data
# The real differences are in how they get applied, but I wanted different type names. 

function MakePeriodicDivergence(mygrid::PeriodicEulGrid)
	# ωx = FFTW.fftfreq(mygrid.Nx,mygrid.Nx/mygrid.L)
	# ωy = FFTW.fftfreq(mygrid.Ny,mygrid.Ny/mygrid.H)
	# compfreqx = im*ωx*(2*pi); #The eigenvalues of a single derivative in x&y
    # compfreqy = im*ωy*(2*pi);
    compfreqx = FirstDerivativeEigenvalues(mygrid.Nx,mygrid.L);
	compfreqy = FirstDerivativeEigenvalues(mygrid.Ny,mygrid.H);


    ΩX = [compfreqx[i] for i=1:mygrid.Nx,j=1:mygrid.Ny]; #put them into arrays the same size as the grid
    ΩY = [compfreqy[j] for i=1:mygrid.Nx,j=1:mygrid.Ny];

	Div = DivergencePeriodicDifferentialOperator(mygrid,ΩX,ΩY);
	return Div
end

function MakePeriodicGradient(mygrid::PeriodicEulGrid)
	# ωx = FFTW.fftfreq(mygrid.Nx,mygrid.Nx/mygrid.L)
	# ωy = FFTW.fftfreq(mygrid.Ny,mygrid.Ny/mygrid.H)
	# compfreqx = im*ωx*(2*pi); #The eigenvalues of a single derivative in x&y
    # compfreqy = im*ωy*(2*pi);
    compfreqx = FirstDerivativeEigenvalues(mygrid.Nx,mygrid.L);
	compfreqy = FirstDerivativeEigenvalues(mygrid.Ny,mygrid.H);

    ΩX = [compfreqx[i] for i=1:mygrid.Nx,j=1:mygrid.Ny]; #put them into arrays the same size as the grid
    ΩY = [compfreqy[j] for i=1:mygrid.Nx,j=1:mygrid.Ny];

	Grad = GradientPeriodicDifferentialOperator(mygrid,ΩX,ΩY);
	return Grad
end

# Now to make ``inverse'' operators. This is only defined for ``simple'' derivatives. 
# Possibly revisit at a later version of the code. 
# Modes whose eigenvalue is (numerically) zero are dropped: the solution gets 0 there.
function InvertOperator(op::SimplePeriodicDifferentialOperator; tol=1e-12)
	cutoff = tol*maximum(abs.(op.Eigenvalues))
	dropped = abs.(op.Eigenvalues) .<= cutoff
	inv = [d ? zero(ComplexF64) : 1/λ for (λ,d) in zip(op.Eigenvalues,dropped)]
	return InverseSimplePeriodicDifferentialOperator(op.grid, inv, dropped)
end


# Making coposition operator div(grad(.)), this is a sort of Laplacian, but 
# built from the zeroed first-derivative eigenvalues (not the true Laplacian).
# Zero at k=0 and, on even grids, on the Nyquist lines; InvertOperator drops those modes.
# Going to use the previously defined operator algebra to do this. 
function MakePeriodicDivGrad(D::DivergencePeriodicDifferentialOperator,G::GradientPeriodicDifferentialOperator)
	return D*G;
end

#Can be made directly from the grid, if D & G don't already exist. 
function MakePeriodicDivGrad(mygrid::PeriodicEulGrid)
	D = MakePeriodicDivergence(mygrid);
	G = MakePeriodicGradient(mygrid);
	return D*G;
end



#Helper function to get the spectral eigenvalues of a 1D derivative 
function FirstDerivativeEigenvalues(N::Integer, L::Real)
	k = im*collect(FFTW.fftfreq(N, N/L))*(2*pi);
	if iseven(N); k[N÷2 + 1] = 0; end
	return k
end

###############################################
#   This is how we apply/invert operators     #
###############################################


#####Simple ones first, there are a few cases#####
#####First we do scalar data, either struct or simple arrays#####
#This applys a simple differential operator to a scalar function on a periodic grid
function ApplySimpleOperator(mydata::ScalarGridData,myoperator::SimplePeriodicDifferentialOperator)
	#`===` is an identity check: data and operator must share the very same grid object,
	#not merely an equal one (see the note on PeriodicEulGrid in EulGrid.jl).
	if ~(mydata.grid === myoperator.grid)
		throw(ArgumentError("Data and operator grids do not match"));
	end
	oldhat = FFTW.fft(mydata.U);
	newhat = myoperator.Eigenvalues.*oldhat;
	new = real(FFTW.ifft(newhat));
	result = ScalarGridData(new,myoperator.grid);
	return result
end

#Same, but no need for a struct. Will apply to an array of data
function ApplySimpleOperator(mydata::Matrix{Float64},myoperator::SimplePeriodicDifferentialOperator)
	if ~(size(mydata) == (myoperator.grid.Nx,myoperator.grid.Ny))
		throw(ArgumentError("Size of input does not match grid size of operator"));
	end
	oldhat = FFTW.fft(mydata);
	newhat = myoperator.Eigenvalues.*oldhat;
	new = real(FFTW.ifft(newhat));
	return new
end


####Now the same as above, but for vector structs or pairs of arrays####

#Applying the operator to a vector grid data struct
function ApplySimpleOperator(mydata::VectorGridData,myoperator::SimplePeriodicDifferentialOperator)
	#Grid identity check (`===`), see ApplySimpleOperator for scalar data.
	if ~(mydata.grid === myoperator.grid)
		throw(ArgumentError("Data and operator grids do not match"));
	end
	oldUhat = FFTW.fft(mydata.U);
	oldVhat = FFTW.fft(mydata.V);
	newUhat = myoperator.Eigenvalues.*oldUhat;
	newVhat = myoperator.Eigenvalues.*oldVhat;
	newU = real(FFTW.ifft(newUhat));
	newV = real(FFTW.ifft(newVhat));
	result = VectorGridData(newU,newV,myoperator.grid);
	return result
end

#Applying it to two arrays
function ApplySimpleOperator(mydataU::Matrix{Float64},mydataV::Matrix{Float64},myoperator::SimplePeriodicDifferentialOperator)
	if ~(size(mydataU) == (myoperator.grid.Nx,myoperator.grid.Ny))
		throw(ArgumentError("Size of first input does not match grid size"));
	elseif ~(size(mydataV) == (myoperator.grid.Nx,myoperator.grid.Ny))
		throw(ArgumentError("Size of second input does not match grid size"));
	end
	oldUhat = FFTW.fft(mydataU);
	oldVhat = FFTW.fft(mydataV);
	newUhat = myoperator.Eigenvalues.*oldUhat;
	newVhat = myoperator.Eigenvalues.*oldVhat;
	newU = real(FFTW.ifft(newUhat));
	newV = real(FFTW.ifft(newVhat));
	return newU,newV
end



#####Now the divergence-like opeators#####
#####First we do vector structs (return scalar struct)#####
function ApplyDivergenceOperator(mydata::VectorGridData,myoperator::DivergencePeriodicDifferentialOperator)
	#Grid identity check (`===`), see ApplySimpleOperator for scalar data.
	if ~(mydata.grid === myoperator.grid)
		throw(ArgumentError("Data and operator grids do not match"));
	end
	oldUhat = FFTW.fft(mydata.U);
	oldVhat = FFTW.fft(mydata.V);
	newUhat = myoperator.EigenvaluesU.*oldUhat;
	newVhat = myoperator.EigenvaluesV.*oldVhat;
	new = real(FFTW.ifft(newUhat .+ newVhat));
	result = ScalarGridData(new,myoperator.grid);
	return result
end

#Now apply to vector data stored as two arrays
function ApplyDivergenceOperator(mydataU::Matrix{Float64},mydataV::Matrix{Float64},myoperator::DivergencePeriodicDifferentialOperator)
	if ~(size(mydataU) == (myoperator.grid.Nx,myoperator.grid.Ny))
		throw(ArgumentError("Size of first input does not match grid size"));
	elseif ~(size(mydataV) == (myoperator.grid.Nx,myoperator.grid.Ny))
		throw(ArgumentError("Size of second input does not match grid size"));
	end
	oldUhat = FFTW.fft(mydataU);
	oldVhat = FFTW.fft(mydataV);
	newUhat = myoperator.EigenvaluesU.*oldUhat;
	newVhat = myoperator.EigenvaluesV.*oldVhat;
	new = real(FFTW.ifft(newUhat .+ newVhat));
	return new
end



#####Next, the gradient-like opeators#####
#####First we do scalar structs (return vector struct)#####
function ApplyGradientOperator(mydata::ScalarGridData,myoperator::GradientPeriodicDifferentialOperator)
	#Grid identity check (`===`), see ApplySimpleOperator for scalar data.
	if ~(mydata.grid === myoperator.grid)
		throw(ArgumentError("Data and operator grids do not match"));
	end
	oldUhat = FFTW.fft(mydata.U);
	newUhat = myoperator.EigenvaluesU.*oldUhat;
	newVhat = myoperator.EigenvaluesV.*oldUhat;
	newU = real(FFTW.ifft(newUhat));
	newV = real(FFTW.ifft(newVhat));
	result = VectorGridData(newU,newV,myoperator.grid);
	return result
end

#Now apply to vector data stored as two arrays
function ApplyGradientOperator(mydata::Matrix{Float64},myoperator::GradientPeriodicDifferentialOperator)
	if ~(size(mydata) == (myoperator.grid.Nx,myoperator.grid.Ny))
		throw(ArgumentError("Size of first input does not match grid size"));
	end
	oldhat = FFTW.fft(mydata);
	newUhat = myoperator.EigenvaluesU.*oldhat;
	newVhat = myoperator.EigenvaluesV.*oldhat;
	newU = real(FFTW.ifft(newUhat));
	newV = real(FFTW.ifft(newVhat));
	return newU, newV
end




##### Inverting Operators simply mirrors apply, but for Inverse Operators ####


#####There are a few cases#####
#####First we do scalar data, either struct or simple arrays#####
#This inverts a simple differential operator to a scalar function on a periodic grid
function ApplyInverseOperator(mydata::ScalarGridData,myoperator::InverseSimplePeriodicDifferentialOperator)
	#`===` is an identity check: data and operator must share the very same grid object,
	#not merely an equal one (see the note on PeriodicEulGrid in EulGrid.jl).
	if ~(mydata.grid === myoperator.grid)
		throw(ArgumentError("Data and operator grids do not match"));
	end
	oldhat = FFTW.fft(mydata.U);
	CheckDroppedModes(oldhat, myoperator);
	newhat = myoperator.InvEigenvalues.*oldhat;
	new = real(FFTW.ifft(newhat));
	result = ScalarGridData(new,myoperator.grid);
	return result
end

#Same, but no need for a struct. Will invert on an array of data
function ApplyInverseOperator(mydata::Matrix{Float64},myoperator::InverseSimplePeriodicDifferentialOperator)
	if ~(size(mydata) == (myoperator.grid.Nx,myoperator.grid.Ny))
		throw(ArgumentError("Size of input does not match grid size of operator"));
	end
	oldhat = FFTW.fft(mydata);
	CheckDroppedModes(oldhat, myoperator);
	newhat = myoperator.InvEigenvalues.*oldhat;
	new = real(FFTW.ifft(newhat));
	return new
end


####Now the same as above, but for vector structs or pairs of arrays####

#Inverting the operator on a vector grid data struct
function ApplyInverseOperator(mydata::VectorGridData,myoperator::InverseSimplePeriodicDifferentialOperator)
	#Grid identity check (`===`), see ApplySimpleOperator for scalar data.
	if ~(mydata.grid === myoperator.grid)
		throw(ArgumentError("Data and operator grids do not match"));
	end
	oldUhat = FFTW.fft(mydata.U);
	oldVhat = FFTW.fft(mydata.V);
	CheckDroppedModes(oldUhat, myoperator);
	CheckDroppedModes(oldVhat, myoperator);
	newUhat = myoperator.InvEigenvalues.*oldUhat;
	newVhat = myoperator.InvEigenvalues.*oldVhat;
	newU = real(FFTW.ifft(newUhat));
	newV = real(FFTW.ifft(newVhat));
	result = VectorGridData(newU,newV,myoperator.grid);
	return result
end

#Inverting it on two arrays
function ApplyInverseOperator(mydataU::Matrix{Float64},mydataV::Matrix{Float64},myoperator::InverseSimplePeriodicDifferentialOperator)
	if ~(size(mydataU) == (myoperator.grid.Nx,myoperator.grid.Ny))
		throw(ArgumentError("Size of first input does not match grid size"));
	elseif ~(size(mydataV) == (myoperator.grid.Nx,myoperator.grid.Ny))
		throw(ArgumentError("Size of second input does not match grid size"));
	end
	oldUhat = FFTW.fft(mydataU);
	oldVhat = FFTW.fft(mydataV);
	CheckDroppedModes(oldUhat, myoperator);
	CheckDroppedModes(oldVhat, myoperator);
	newUhat = myoperator.InvEigenvalues.*oldUhat;
	newVhat = myoperator.InvEigenvalues.*oldVhat;
	newU = real(FFTW.ifft(newUhat));
	newV = real(FFTW.ifft(newVhat));
	return newU,newV
end



#Helper function to test that we're not throwing away data 
# Warn if the transformed right-hand side has content in modes the operator cannot invert.
# For Poisson that is the mean (and any Nyquist lines for D*G); that content is discarded.
function CheckDroppedModes(fhat, myoperator::InverseSimplePeriodicDifferentialOperator; tol=1e-8)
	total = sqrt(sum(abs2, fhat))
	lost  = sqrt(sum(abs2, fhat[myoperator.dropped]))
	if lost > tol*total
		@warn "Right-hand side has content in modes the operator cannot invert (for Poisson: a nonzero mean). That content is discarded." relative_size=lost/total maxlog=10
	end
end




# Helper Function to return a boolean determining if two operators act on the same grid
function CheckSameGrid(A::AbstractDifferentialOperator,B::AbstractDifferentialOperator,what::String)
	if ~(A.grid === B.grid)
		throw(ArgumentError("Operators to be $what are on different grids"));
	end
end