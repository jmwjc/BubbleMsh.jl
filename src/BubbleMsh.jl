module BubbleMsh

"""
Bubble mesh generator implemented by `Julia`
"""

import Triangle: basic_triangulation
export bubblemsh

"""
bubblemsh(filename::String,xc::Vector{Float64},d::Vector{Float64},n::Int,r::Float64,h::Float64)
Input:
    filename: name of `msh` file 
    xc: the coordiantes of centric point used for generating
    d: the max distance between centric point and generated point in each directions
    n: number of generated points
    r: characteristic distance between each points
    h: length of time steps in Runge-Kutta method
Output:
    a `msh` file named by original filename with a suffix of `bubble`
"""
function bubblemsh(filename::String,xc::Vector{Float64},d::Vector{Float64},n::Int,r::Float64,h::Float64)
    fi = open(filename,"r")
    filename = replace(filename,r"\.msh"=>s"_bubble.msh")
    fo = open(filename,"w")
    k = 0
    phytag = 0
    x = zeros(n,3)
    for line in eachline(fi)
        if line == "\$PhysicalNames"
            write(fo,line*"\n")
            m = parse(Int,readline(fi))
            write(fo,"$(m+1)\n")
            for i in 1:m
                write(fo,readline(fi)*"\n")
            end
            phytag = m+1
            write(fo,"$phytag 2 \"Ω\"\n")
        elseif line == "\$Nodes"
            write(fo,line*"\n")
            k = parse(Int,readline(fi))
            write(fo,"$(n+k)\n")
            xᵇ = zeros(k,3)
            for i in 1:k
                line = readline(fi)
                write(fo,line*"\n")
                t_,x_,y_,z_ = split(line," ")
                tag = parse(Int,t_)
                xᵇ[i,1] = parse(Float64,x_)
                xᵇ[i,2] = parse(Float64,y_)
                xᵇ[i,3] = parse(Float64,z_)
            end
            x .= initnodes(xc,d,n)
            bubble!(x,xᵇ,r,h)
            for i in 1:n
                write(fo,"$(k+i) $(x[i,1]) $(x[i,2]) $(x[i,3])\n")
            end
            x = vcat(xᵇ,x)
        elseif line == "\$Elements"
            write(fo,line*"\n")
            m = parse(Int,readline(fi))
            write(fo,"$(2*(n+k)-2-k+m)\n")
            for i in 1:m
                write(fo,readline(fi)*"\n")
            end
            triangle = basic_triangulation(x,collect(1:(n+k)))
            for (i,ids) in enumerate(triangle)
                write(fo,"$(m+i) 2 2 $phytag 1 $(ids[1]) $(ids[2]) $(ids[3])\n")
            end
        else
            write(fo,line*"\n")
        end
    end
    close(fi)
    close(fo)
end

"""
bubble
The interface force between two neighbor bubbles is defined by Dinh, V. Q., & Marechal, Y. (2017). GPU-based parallelization for bubble mesh generation. COMPEL - The International Journal for Computation and Mathematics in Electrical and Electronic Engineering, 36(4), 1184–1197. https://doi.org/10.1108/COMPEL-11-2016-0476
"""
function bubble!(x::Matrix{Float64},x_::Matrix{Float64},r::Float64,h::Float64)
    n = size(x_,1)
    nₚ = size(x,1)
    k = 1.47/2/r
    c = 1.4*k^0.5
    iter = 0
    while iter < 10000
        iter += 1
        p = zeros(nₚ,3)
        for i in 1:nₚ
            xᵢ = x[i,1]
            yᵢ = x[i,2]
            zᵢ = x[i,3]
            for j in 1:n
                xⱼ = x_[j,1]
                yⱼ = x_[j,2]
                zⱼ = x_[j,3]
                Δx = xᵢ-xⱼ
                Δy = yᵢ-yⱼ
                Δz = zᵢ-zⱼ
                dij = (Δx^2+Δy^2+Δz^2)^0.5
                w = dij/r/2
                Δp = w>1.5 ? 0.0 : (1.0-w^4)*exp(-w^4)
                p[i,1] += Δp*Δx/dij
                p[i,2] += Δp*Δy/dij
                p[i,3] += Δp*Δz/dij
            end
            for j in 1:nₚ
                if i ≠ j
                    xⱼ = x[j,1]
                    yⱼ = x[j,2]
                    zⱼ = x[j,3]
                    Δx = xᵢ-xⱼ
                    Δy = yᵢ-yⱼ
                    Δz = zᵢ-zⱼ
                    dij = ((xᵢ-xⱼ)^2+(yᵢ-yⱼ)^2+(zᵢ-zⱼ)^2)^0.5
                    w = dij/r/2
                    Δp = w>1.5 ? 0.0 : (1.0-w^4)*exp(-w^4)
                    p[i,1] += Δp*Δx/dij
                    p[i,2] += Δp*Δy/dij
                    p[i,3] += Δp*Δz/dij
                end
            end
        end

        rep = p[:,1]'*p[:,1] + p[:,2]'*p[:,2] + p[:,3]'*p[:,3]
        println("iter = $iter")
        println("rep = $rep")
        rep < 1e-10 ? break : nothing

        y = zeros(nₚ,3)
        for i in 1:nₚ
            for j in 1:3
                yᵢ = y[i,j]

                k₁ = h*(p[i,j]-c*yᵢ)
                k₂ = h*(p[i,j]-c*(yᵢ+k₁/2))
                k₃ = h*(p[i,j]-c*(yᵢ+k₂/2))
                k₄ = h*(p[i,j]-c*(yᵢ+k₃))
                y[i,j] = yᵢ + (k₁+2*k₂+2*k₃+k₄)/6

                yᵢ = y[i,j]
                k₁ = h*yᵢ
                k₂ = h*(yᵢ+k₁/2)
                k₃ = h*(yᵢ+k₂/2)
                k₄ = h*(yᵢ+k₃)
                x[i,j] += (k₁+2*k₂+2*k₃+k₄)/6
            end
        end
    end
end

function initnodes(x::Vector{Float64},d::Vector{Float64},n::Int)
    x_ = zeros(n,3)
    for i in 1:n
        x_[i,1] = x[1] + (2*rand()-1.0)*d[1]
        x_[i,2] = x[2] + (2*rand()-1.0)*d[2]
        x_[i,3] = x[3] + (2*rand()-1.0)*d[3]
    end
    return x_
end

end
