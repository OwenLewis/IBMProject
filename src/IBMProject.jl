module IBMProject

# Custom types
include("EulGrid.jl")
# include("ImmersedBnd.jl")
# include("ImmersedMesh.jl")
# include("PDETypes.jl")

# # Mesh input/generation
# include("MeshGeneration.jl")
# include("GMSHreader.jl")
# include("MeshTransform.jl")

# # PostProcessing
# include("vtkExport.jl")
# include("ErrorAnalysis.jl")


# # EulGrid.jl
export PeriodicEulGrid,
       ScalarGridData,
       VectorGridData,
       ApplySingleOperator,
       InvertSingleOperator
       # BrinkmanParam,
       # BrinkmanMPParam,
       # FluidVarParam

# # ImmersedBnd.jl
# export Problem,
#        Dirichlet,
#        Neumann,
#        Robin,
#        Forcing

# # ImmersedMesh.jl
# export squareMesh,
#        squareMeshFluid,
#        Mesh,
#        FluidMesh,
#        axisymMesh,
#        axisymFluid

# # PDETypes.jl
# export pointTransform,
#        onSegment,
#        doIntersect,
#        isInsideDomain,
#        meshTransform

# # ErrorAnalysis.jl
# export DomainNorm,
#        hCalc
end # module

