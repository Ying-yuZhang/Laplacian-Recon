include("this.jl")
using LinearAlgebra, DelimitedFiles, MAT
input = matread("this_input.mat")
X = input["XX"]
Y = input["YY"]
# 超图推断参数
ooi = [1,2]
dmax = 2

# 使用 THIS 推断
Ainf, coeff, relerr = this(X, Y, ooi, dmax,1e-2,1e-2)

mat2=Ainf[2]
this_EdgeList = mat2[:, 1:3]

mat3=Ainf[3]
this_TriangleList = mat3[:, 1:4]
# -------------------------
# 保存结果为 .mat 文件
# -------------------------
matwrite("this_result.mat", Dict(
    "this_EdgeList" => this_EdgeList,
    "this_TriangleList" => this_TriangleList
))