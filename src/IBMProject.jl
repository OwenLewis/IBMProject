module IBMProject

# Custom types
include("EulGrid.jl")
include("ImmersedBnd.jl")
include("DiffOpers.jl")
include("ImmersedNetwork.jl")

# # Network Functions
# include("NetworkOpers.jl")


# # PostProcessing
# include("vtkExport.jl")
# include("ErrorAnalysis.jl")


# # EulGrid.jl
export PeriodicEulGrid,
       ScalarGridData,
       VectorGridData,
       GridIntegral


# ImmersedBnd.jl
export PeriodicLagBnd,
       ScalarBndData,
       ScalarBndSpread,
       ScalarBndInterp,
       BndIntegral

# DiffOpers.jl
export SimplePeriodicDifferentialOperator,
       GradientPeriodicDifferentialOperator,
       DivergencePeriodicDifferentialOperator,
       ApplySimpleOperator,
       ApplyDivergenceOperator,
       ApplyGradientOperator


#ImmersedNetwork.jl
export MeshVertex,
       MeshEdge,
       MeshFace, 
       LagMesh

end # module

