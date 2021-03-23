# ------------------------------------
# set initial conditions
# ------------------------------------
function set_initQbase(Qbase, cellxmax, cellymax, cellzmax, restart_file, init_rho, init_u, init_v, init_w, init_p, init_T,
                    specific_heat_ratio, out_file_front, out_ext, out_dir, restartnum, Rd, nval, icell)
    
    # 角セルに値を代入するため
    Qbase = setup_init_value(Qbase, cellxmax, cellymax, cellzmax, init_rho, init_u, init_v, init_w, init_p)

    restart_check = 0
    try Qbase = setup_restart_value(Qbase, cellxmax, cellymax, cellzmax, out_dir, restart_file, nval, icell)
        println("Restart "*restart_file)
        restart_check = 2
    catch 
        restart_check = 1
    end

    if restart_check == 1
        Qbase = setup_init_value(Qbase, cellxmax, cellymax, cellzmax, init_rho, init_u, init_v, init_w, init_p)
        println("Start Initial condition")
        restartnum = 0
        output_result(0, Qbase, cellxmax, cellymax, cellzmax, specific_heat_ratio, out_file_front, out_ext, out_dir, Rd, nval, icell)
    end

    return Qbase, restartnum
end

function setup_init_value(Qbase, cellxmax, cellymax, cellzmax, init_rho, init_u, init_v, init_w, init_p)
    for k in 1:cellzmax
        for j in 1:cellymax
            for i in 1:cellxmax
                Qbase[i,j,k,1] = init_rho
                Qbase[i,j,k,2] = init_u
                Qbase[i,j,k,3] = init_v
                Qbase[i,j,k,4] = init_w
                Qbase[i,j,k,5] = init_p
            end
        end
    end
    return Qbase
end

function setup_restart_value(Qbase, cellxmax, cellymax, cellzmax, out_dir, restart_file, nval, icell)
    skipnum = 1
    fff = []
    open("result/"*restart_file, "r") do f
        fff = read(f, String)
    end 
    fff = split(fff,"\n", keepempty = false)
    
    for i in 1+skipnum:length(fff)
        fff[i] = replace(fff[i]," \r" => "")
    end
    
    ite = 1
    for i in 1+icell:cellxmax-icell
        for j in 1+icell:cellymax-icell
            for k in 1+icell:cellzmax-icell
                temp = split(fff[ite+skipnum]," ")
                for l in 1:nval
                    Qbase[i,j,k,l] = parse(Float64,temp[l]) 
                end
                ite = ite+1
            end
        end
    end
    return Qbase
end

# ------------------------------------
# find out if results were diverge
# ------------------------------------
function check_divrege(Qbase, cellxmax, cellymax, cellzmax, Rd, fwrite, icell)
    ite = 0
    for k in 1+icell:cellzmax-icell
        for j in 1+icell:cellymax-icell
            for i in 1+icell:cellxmax-icell        
                T = Qbase[i,j,k,5]/(Qbase[i,j,k,1]*Rd)
                if T < 0 || isequal(T, NaN) == true
                    open( fwrite, "a" ) do f
                        ai  = @sprintf("%4.0f", i)
                        aj  = @sprintf("%4.0f", j)
                        ak  = @sprintf("%4.0f", k)
                        rho = @sprintf("%4.0e", Qbase[i,j,k,1])
                        p   = @sprintf("%4.0e", Qbase[i,j,k,5])

                        write(f, "\n")
                        write(f, " diverge ")
                        write(f, "\n")
                        write(f, " i = "*ai)
                        write(f, "\n")
                        write(f, " j = "*aj)
                        write(f, "\n")
                        write(f, " k = "*ak)
                        write(f, "\n")
                        write(f, " rho = "*rho)
                        write(f, "\n")
                        write(f, " p = "*p)
                        write(f, "\n")

                    end
                    ite = 1
                end
            end
        end
    end

    if ite == 1
        println("\n")
        println("\n")
        println(" T<0 ")
        println(" diverge ")
        println("\n")
        throw(UndefVarError(:x))
    end
end

# ------------------------------------
# quicksort
# ------------------------------------
function quicksort(list, first, last, num_cell_data, standard_num)
    # list = zezros(m,n)
    x = list[Int(floor((first+last)/2)),standard_num]
    i = first
    j = last
    
    while true
        while list[i,standard_num] < x
            i = i+1
        end
        while x < list[j,standard_num]
            j = j-1
        end
        
        if (i >= j) 
            break
        else
            for k in 1:num_cell_data
                t = list[i,k]
                list[i,k] = list[j,k]
                list[j,k] = t
            end
            i = i+1
            j = j-1
        end
    end
    
    if (first < i - 1) 
        list = quicksort(list, first, i - 1,num_cell_data,standard_num)
    end
    if (j + 1 < last) 
        list = quicksort(list, j + 1, last,num_cell_data,standard_num)
    end
    return list
end