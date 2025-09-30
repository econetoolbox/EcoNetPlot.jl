using CSV, DataFrames, Serialization

df = CSV.read("data/chilean_TI.txt", DataFrame)
select!(df, Not(:Column1, :Column2))
A_trophic = Array(transpose(Matrix(df)))
serialize("data/A_trophic.bin", A_trophic)


df_pos = CSV.read("data/chilean_NTIpos.txt", DataFrame)
select!(df_pos, Not(:Column1, :Column2))
A_pos = Array(transpose(Matrix(df_pos)))
serialize("data/A_positive.bin", A_pos)

df_neg = CSV.read("data/chilean_NTIneg.txt", DataFrame)
select!(df_neg, Not(:Column1, :Column2))
A_neg = Array(transpose(Matrix(df_neg)))
serialize("data/A_negative.bin", A_neg)

