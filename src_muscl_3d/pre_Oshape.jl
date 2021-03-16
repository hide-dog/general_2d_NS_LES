using Printf

function main()

    infile = "xy_hayabusa"
    znum = 5
    dz = 1e-3

    xnum, ynum, x, y = read_wing(infile)

    outdir = "grid"
    make_dir(outdir)
    result = "result"
    make_dir(result)
    result = "post_result"
    make_dir(result)
    nodes, xnum_max, ynum_max, znum_max = mk_gird(xnum, ynum, znum, x, y, dz, outdir)
    vecA(nodes, xnum_max, ynum_max, znum_max, outdir)
end

function read_wing(infile)
    
    fff=[]
    open(infile, "r") do f  # 全行格納
        fff = read(f,String)
    end 
    fff = split(fff,"\n",keepempty=false)   # 改行分割(\n)
    
    for i in 1:length(fff)
        fff[i] = replace(fff[i]," \r" => "")  # 改行削除(\r)
    end

    temp   = split(fff[1]," ")
    filter!(e->e≠"",temp)

    xnum = Int(parse(Float64,temp[1]))
    ynum = Int(parse(Float64,temp[2]))
    
    x = zeros(xnum, ynum)
    y = zeros(xnum, ynum)
    skip = 1
    fnum = length(fff)-skip
    for i in 1:fnum
        temp = split(fff[i+skip]," ")
        filter!(e->e≠"",temp)
        
        tx = parse(Int64,temp[1])
        ty = parse(Int64,temp[2])
        x[tx,ty] = parse(Float64,temp[3])
        y[tx,ty] = parse(Float64,temp[4])
    end

    return xnum, ynum, x, y 
end

function mk_gird(xnum, ynum, znum, x, y, dz, outdir)
    """
    nodes[i,j,k,l]
    i   : x,r方向の番号
    j   : y,theta方向の番号
    k   : z方向の番号
    l=1 : 点のx座標
    l=2 : 点のy座標

    nodes[1,:]      : x方向境界
    nodes[:,1]      : y方向境界
    nodes[xnum+2,:] : x方向境界
    nodes[:,ynum+2] : y方向境界
    """
    icell = 2
    xnum_max = xnum + 2*icell
    ynum_max = ynum + 2*icell
    znum_max = znum + 2*icell
    nodes = zeros(xnum_max, ynum_max, znum_max, 3)

    
    for j in 1+icell:ynum_max-icell
        for i in 1+icell:xnum_max-icell
            nodes[i,j,1,1] = x[i-icell,j-icell]
            nodes[i,j,1,2] = y[i-icell,j-icell]
        end
    end

    #=
    仮想セルの作成：現在あるセル境界線を延長して作成
    そのまま延長した際にセルが交差するのを防ぐため、仮想セルの大きさを100^-1倍(n)にする.
    ベクトルで考えれば下記のような計算になる．たぶん
    また、[1,1]等の角にはダミー値が入っている.
    =#
    n = 1.0
    for j in 1:ynum_max
        for l in 1:2
            nodes[2,j,1,l]        =  (1.0+n)*nodes[3,j,1,l] - n*nodes[4,j,1,l]
            nodes[1,j,1,l]        =  (1.0+n)*nodes[2,j,1,l] - n*nodes[3,j,1,l]
            nodes[xnum_max-1,j,1,l] =  (1.0+n)*nodes[xnum_max-2,j,1,l] - n*nodes[xnum_max-3,j,1,l]
            nodes[xnum_max,j,1,l] =  (1.0+n)*nodes[xnum_max-1,j,1,l] - n*nodes[xnum_max-2,j,1,l]
        end
    end
    for i in 1:xnum_max
        for l in 1:2
            nodes[i,2,1,l]        =  (1.0+n)*nodes[i,3,1,l] - n*nodes[i,4,1,l]
            nodes[i,1,1,l]        =  (1.0+n)*nodes[i,2,1,l] - n*nodes[i,3,1,l]
            nodes[i,ynum_max-1,1,l] =  (1.0+n)*nodes[i,ynum_max-2,1,l] - n*nodes[i,ynum_max-3,1,l]
            nodes[i,ynum_max,1,l] =  (1.0+n)*nodes[i,ynum_max-1,1,l] - n*nodes[i,ynum_max-2,1,l]
        end
    end

    for k in 1:znum_max
        for j in 1:ynum_max
            for i in 1:xnum_max
                nodes[i,j,k,1] = nodes[i,j,1,1]
                nodes[i,j,k,2] = nodes[i,j,1,2]
                nodes[i,j,k,3] = dz*(k-1)
            end
        end
    end

    fff = outdir*"/nodes"
    open(fff,"w") do f
        write(f,"nodes: xnum, ynum , x, y\n")
        for i in 1:xnum_max
            for j in 1:ynum_max
                for k in 1:znum_max
                    x = @sprintf("%8.8e", nodes[i,j,k,1])
                    y = @sprintf("%8.8e", nodes[i,j,k,2])
                    z = @sprintf("%8.8e", nodes[i,j,k,3])
                    write(f,string(i)*" "*string(j)*" "*string(k)*" "*x*" "*y*" "*z*"\n")
                end
            end
        end
    end
    println("write "*fff)

    # nodes_forvtk
    nodes_num = zeros(Int, xnum_max, ynum_max, znum_max)
    fff=outdir*"/nodes_forvtk"
    open(fff,"w") do f
        write(f,"nodes: xnum, ynum, znum, x, y, z\n")
        k=1
        for i in 1+icell:xnum_max-icell
            for j in 1+icell:ynum_max-icell
                x = @sprintf("%8.8e", nodes[i,j,1])
                y = @sprintf("%8.8e", nodes[i,j,2])
                z = @sprintf("%8.8e", nodes[i,j,3])
                write(f,string(k)*" "*x*" "*y*" "*z*"\n")
                nodes_num[i,j] = k
                k = k+1
            end
        end
    end
    println("write "*fff)

    fff=outdir*"/nodesnum"
    open(fff,"w") do f
        write(f,"nodesnum: xnum_max, ynum_max, znum_max\n")
        write(f,string(xnum_max) * " "*string(ynum_max) * " "*string(znum_max)*"\n")
    end

    ### element ###
    fff=outdir*"/element_forvtk"
    open(fff,"w") do f
        write(f,"elements:cell_xnum, lup,rup,ldown,rdown, lup,rup,ldown,rdown \n")
        k=1
        for i in 1+icell:xnum_max-icell-1
            for j in 1+icell:ynum_max-icell-1
                for k in 1+icell:znum_max-icell-1
                    d1 = @sprintf("%1.0f", nodes_num[  i,  j,  k])
                    d2 = @sprintf("%1.0f", nodes_num[  i,j+1,  k])
                    d3 = @sprintf("%1.0f", nodes_num[i+1,j+1,  k])
                    d4 = @sprintf("%1.0f", nodes_num[i+1,  j,  k])
                    d5 = @sprintf("%1.0f", nodes_num[  i,  j,k+1])
                    d6 = @sprintf("%1.0f", nodes_num[  i,j+1,k+1])
                    d7 = @sprintf("%1.0f", nodes_num[i+1,j+1,k+1])
                    d8 = @sprintf("%1.0f", nodes_num[i+1,  j,k+1])
                    write(f, string(k)*" "*d1*" "*d2*" "*d3*" "*d4)
                    write(f, " "*d5*" "*d6*" "*d7*" "*d8*"\n")
                    k = k+1
                end
            end
        end
    end
    println("write "*fff)

    return  nodes, xnum_max, ynum_max, znum_max
end

function vecA(nodes,xnum_max,ynum_max,znum_max,outdir)

    vecAx = zeros(xnum_max, ynum_max-1, znum_max-1, 3)
    for i in 1:xnum_max
        for j in 1:ynum_max-1
            for k in 1:znum_max-1
                # 3dim Ax=(, , 0)
                x = nodes[i,j+1,k,2] - nodes[i,j,k,2]
                y = nodes[i,j,k,1] - nodes[i,j+1,k,1]
                z = nodes[i,j,k,3] - nodes[i,j+1,k,3]
                vecAx[i,j,k,1] = x
                vecAx[i,j,k,2] = y
                vecAx[i,j,k,3] = z
            end
        end
    end

    vecAy = zeros(xnum_max-1, ynum_max, znum_max-1, 3)
    for i in 1:xnum_max-1
        for j in 1:ynum_max
            for k in 1:znum_max-1
            # 3dim Ay=( , , 0)
            x = nodes[i,j,k,2]-nodes[i+1,j,k,2]
            y = nodes[i+1,j,k,1]-nodes[i,j,k,1]
            z = nodes[i+1,j,k,1]-nodes[i,j,k,1]
            vecAy[i,j,k,1] = x
            vecAy[i,j,k,2] = y
            vecAy[i,j,k,3] = z
        end
    end

    vecAz = zeros(xnum_max-1, ynum_max-1, znum_max, 3)
    for i in 1:xnum_max-1
        for j in 1:ynum_max-1
            for k in 1:znum_max
            # 3dim Az=( , , 0)
            x = nodes[i,j,k,2]-nodes[i+1,j,k,2]
            y = nodes[i+1,j,k,1]-nodes[i,j,k,1]
            z = nodes[i+1,j,k,1]-nodes[i,j,k,1]
            vecAz[i,j,k,1] = x
            vecAz[i,j,k,2] = y
            vecAz[i,j,k,3] = z
        end
    end

    fff = outdir*"/vecAx"
    open(fff,"w") do f
        write(f,"vecAx: xnum, ynum, znum, x vec, y vec, z vec\n")
        for i in 1:xnum_max
            for j in 1:ynum_max-1
                for k in 1:znum_max-1
                    x = @sprintf("%8.8f", vecAx[i,j,k,1])
                    y = @sprintf("%8.8f", vecAx[i,j,k,2])
                    z = @sprintf("%8.8f", vecAx[i,j,k,3])
                    write(f, string(i)*" "*string(j)*" "*string(k)*" "*x*" "*y*" "*z*"\n")
                end
            end
        end
    end
    println("write "*fff)

    fff = outdir*"/vecAy"
    open(fff,"w") do f
        write(f,"vecAy: xnum, ynum, znum, x vec, y vec, z vec\n")
        for i in 1:xnum_max-1
            for j in 1:ynum_max
                for k in 1:znum_max-1
                    x = @sprintf("%8.8f", vecAy[i,j,k,1])
                    y = @sprintf("%8.8f", vecAy[i,j,k,2])
                    z = @sprintf("%8.8f", vecAy[i,j,k,3])
                    write(f, string(i)*" "*string(j)*" "*string(k)*" "*x*" "*y*" "*z*"\n")
                end
            end
        end
    end
    println("write "*fff)

    fff = outdir*"/vecAz"
    open(fff,"w") do f
        write(f,"vecAz: xnum, ynum, znum, x vec, y vec, z vec\n")
        for i in 1:xnum_max-1
            for j in 1:ynum_max-1
                for k in 1:znum_max
                    x = @sprintf("%8.8f", vecAz[i,j,k,1])
                    y = @sprintf("%8.8f", vecAz[i,j,k,2])
                    z = @sprintf("%8.8f", vecAz[i,j,k,3])
                    write(f, string(i)*" "*string(j)*" "*string(k)*" "*x*" "*y*" "*z*"\n")
                end
            end
        end
    end
    println("write "*fff)
end

function make_dir(outdir)
    k = 0
    try rm(outdir,recursive=true)
    catch
        mkdir(outdir)
        k = 1
    end

    if k == 0
        mkdir(outdir)
    end
end

# ---------------------------------
main() 