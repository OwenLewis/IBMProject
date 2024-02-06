module IBMProject

# Custom types
include("EulGrid.jl")
include("ImmersedBnd.jl")
include("DiffOpers.jl")
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
       PeriodicDifferentialOperator,
       ApplySingleOperator,
       InvertSingleOperator


# ImmersedBnd.jl
export PeriodicLagBnd,
       ScalarBndData,
       ScalarBndSpread,
       ScalarBndInterp,
       BndIntegral


end # module

