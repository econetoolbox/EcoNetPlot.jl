using CSV, DataFrames

df = CSV.read("data/chilean_TI.txt", DataFrame)
select!(df, Not(:Column1, :Column2))
A_trophic = Array(transpose(Matrix(df)))

df_pos = CSV.read("data/chilean_NTIpos.txt", DataFrame)
select!(df_pos, Not(:Column1, :Column2))
A_pos = Array(transpose(Matrix(df_pos)))

df_neg = CSV.read("data/chilean_NTIneg.txt", DataFrame)
select!(df_neg, Not(:Column1, :Column2))
A_neg = Array(transpose(Matrix(df_neg)))

using EcologicalNetworksDynamics, EcoNetPlot

fw = Foodweb(A_trophic)
plot_network(fw; A_facilitation = A_pos, A_competition = A_neg)
