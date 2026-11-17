module IBMProject

# Custom types
include("LagGeometry.jl")
include("EulGrid.jl")
include("ImmersedBnd.jl")
include("DiffOpers.jl")
include("ImmersedNetwork.jl")

# # Immersed Bnd & Network Functions
include("NetworkOpers.jl")
include("BndOpers.jl")


# # PostProcessing
# include("vtkExport.jl")
# include("ErrorAnalysis.jl")


# EulGrid.jl
export PeriodicEulGrid,
       ScalarGridData,
       VectorGridData,
       GridIntegral


# ImmersedBnd.jl
export PeriodicLagBnd,
       ScalarBndData,
       ScalarBndSpread,
       ScalarBndInterp,
       VectorBndData,
       BndIntegral

# DiffOpers.jl
export SimplePeriodicDifferentialOperator,
       GradientPeriodicDifferentialOperator,
       DivergencePeriodicDifferentialOperator,
       InverseSimplePeriodicDifferentialOperator,
       MakePeriodicGradient,
       MakePeriodicDivergence,
       MakePeriodicLaplacian,
       MakePeriodicHelmholtz,
       InvertOperator,
       ApplySimpleOperator,
       ApplyDivergenceOperator,
       ApplyGradientOperator,
       ApplyInverseOperator


#ImmersedNetwork.jl
export MeshVertex,
       MeshEdge,
       MeshFace, 
       LagMesh

#NetworkOpers.jl
export ReadNetwork,
       MoveNetwork!,
       ElasticForceCalc

# BndOpers.jl
export ReadBnd,
       MoveBnd!

end # module

