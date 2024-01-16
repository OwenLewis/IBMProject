# Code to be loaded into IBMProject

abstract type AbstractGrid end

struct PeriodicGrid <: AbstractGrid
	Lx::Float64
	Ly::Float64
	Ux::Float64
	Uy::Float64
	Nx::Int
	Ny::Int
	dx::Float64
	dy::Float64
	L::Float64
	H::Float64
end


#No need to specify the redundant information
function PeriodicGrid(Lx::Real,Ly::Real,Ux::Real,Uy::Real,Nx::Int,Ny::Int)
	L = Ux - Lx;
	H = Uy - Ly;
	dx = L/Nx;
	dy = H/Ny;
	grid::PeriodicGrid = PeriodicGrid(Lx,Ly,Ux,Uy,Nx,Ny,dx,dy,L,H)
	return grid
end

#If lower bounds are not specified, we assume they're zero
function PeriodicGrid(Ux::Real,Uy::Real,Nx::Int,Ny::Int)
	Lx = 0.0;
	Ly = 0.0;
	L = Ux - Lx;
	H = Uy - Ly;
	dx = L/Nx;
	dy = H/Ny;
	grid::PeriodicGrid = PeriodicGrid(Lx,Ly,Ux,Uy,Nx,Ny,dx,dy,L,H)
	return grid
end

#If upper bounds are not specified, we assume they're unity
function PeriodicGrid(Nx::Int,Ny::Int)
	Lx = 0.0;
	Ly = 0.0;
	Ux = 1.0;
	Uy = 1.0;
	L = Ux - Lx;
	H = Uy - Ly;
	dx = L/Nx;
	dy = H/Ny;
	grid::PeriodicGrid = PeriodicGrid(Lx,Ly,Ux,Uy,Nx,Ny,dx,dy,L,H)
	return grid
end