using Printf

# ------------------------------------
# output results
# ------------------------------------
function output_result(stepnum, Qbase, cellxmax, cellymax, cellzmax, specific_heat_ratio, 
                        out_file_front, out_ext, out_dir, Rd, nval, icell)
    
    stepnum = string(stepnum)
    while length(stepnum) < 6
        stepnum = "0"*stepnum
    end
    
    fff = out_dir*"/"*out_file_front*stepnum*out_ext
    open(fff,"w") do f
        write(f,"result:rho[kg/m^3], u[m/s], v[m/s], w[m/s], p[Pa], T[K]\n")
        for i in 1+icell:cellxmax-icell
            for j in 1+icell:cellymax-icell
                for k in 1+icell:cellzmax-icell
                    for l in 1:nval
                        a = @sprintf("%8.8e", Qbase[i,j,k,l])
                        write(f, a*" ")
                    end
                    T = Qbase[i,j,k,5]/(Qbase[i,j,k,1]*Rd)
                    a = @sprintf("%8.8e", T)
                    write(f, a*"\n")
                end
            end
        end
    end
    println("\nwrite "*fff)
end


# ------------------------------------
# output results ( + yplus) 
# ------------------------------------
function output_result1(stepnum, Qbase, yplus, cellxmax, cellymax, cellzmax, specific_heat_ratio,
                        out_file_front, out_ext, out_dir, Rd, nval, icell)
    
    stepnum = string(stepnum)
    while length(stepnum) < 6
        stepnum = "0"*stepnum
    end
    
    fff = out_dir*"/"*out_file_front*stepnum*out_ext
    open(fff,"w") do f
        write(f,"result:rho[kg/m^3], u[m/s], v[m/s], w[m/s], p[Pa], T[K], yplus\n")
        for i in 1+icell:cellxmax-icell
            for j in 1+icell:cellymax-icell
                for k in 1+icell:cellzmax-icell
                    for l in 1:nval
                        a = @sprintf("%8.8e", Qbase[i,j,k,l])
                        write(f, a*" ")
                    end
                    T = Qbase[i,j,k,5]/(Qbase[i,j,k,1]*Rd)
                    a = @sprintf("%8.8e", T)
                    write(f, a* " ")
                    a = @sprintf("%8.8e", yplus[i,j,k])
                    write(f, a*"\n")
                end
            end
        end
    end
    println("\nwrite "*fff)
end

# ------------------------------------
# output ave results
# ------------------------------------
function output_ave(Qbase_ave, cellxmax, cellymax, cellzmax, out_file_front, out_ext, out_dir, Rd, nval, loop_ite, icell)

    fff = out_dir*"/"*out_file_front*"average"*out_ext
    open(fff,"w") do f
        write(f,"result:rho[kg/m^3], u[m/s], v[m/s], w[m/s], p[Pa], T[K]\n")
        for i in 1+icell:cellxmax-icell
            for j in 1+icell:cellymax-icell
                for k in 1+icell:cellzmax-icell
                    for l in 1:nval
                        a = @sprintf("%8.8e", Qbase_ave[i,j,k,l] / loop_ite)
                        write(f, a*" ")
                    end
                    T = Qbase_ave[i,j,k,5]/(Qbase_ave[i,j,k,1]*Rd) / loop_ite
                    a = @sprintf("%8.8e", T)
                    write(f, a*"\n")
                end
            end
        end
    end
    println("\nwrite "*fff)
end

function reset_write(fwrite)
    open( fwrite, "w" ) do f 
    end
end

function output_physicaltime(fwrite, t, dt)
    a = string(t)
    b = @sprintf("%8.8f", t*dt)
    open(fwrite, "a") do f
        write(f, "\n")
        write(f, "step = "*a)
        write(f, "     ")
        write(f, "time[s] = "*b)
        write(f, "\n")
        write(f, "\n")
    end
end

function output_innertime(fwrite, tau, norm2, nval)
    a = string(tau)
    open(fwrite, "a") do f
        write(f, a)
        write(f, " ")
        for l in 1:nval
            a = @sprintf("%8.8e", norm2[l])
            write(f, a)
            write(f, " ")
        end
        write(f, "\n")
    end
end

function output_fin(fwrite, start_t, end_t, nt, dt, in_nt, cellxmax, cellymax, cellzmax)
    
    etime = end_t - start_t     # 経過時間
    outtime = ["temp"]
    outtime[1] = string(etime)        # 文字列への変換

    et = replace(outtime[1], " milliseconds"=>"")
    et = parse(Float64, et)           # 整数への変換(ミリ秒)

    et_h = floor(et / 10^3 / 60 /60)  # n時間の計算
    et = et - et_h * (10^3 * 60 * 60) # et - n hour
    et_m = floor(et / 10^3 / 60 )     # n分の計算
    et = et - et_m * (10^3 * 60)      # et - n min
    et_s = floor(et / 10^3)           # n秒の計算
    et = et - et_s * (10^3)           # et -n sec

    et_h = @sprintf( "%02d", et_h)    # 文字列への変換
    et_m = @sprintf( "%02d", et_m)    # 文字列への変換
    et_s = @sprintf( "%02d", et_s)    # 文字列への変換
    et = et_h*":"*et_m*":"*et_s       # 結合

    numcell = string(cellxmax * cellymax * cellzmax)
    dt      = string(dt)
    nt      = string(nt)
    in_nt   = string(in_nt)

    open(fwrite, "a") do f
        write(f, "\n")
        write(f, "\n")
        write(f, "---------------------------- \n")
        write(f, "  fin program  ")
        write(f, "\n")
        write(f, "  number of cell            = " * numcell)
        write(f, "\n")
        write(f, "  dt                        = " * dt)
        write(f, "\n")
        write(f, "  number of time iteration  = " * nt)
        write(f, "\n")
        write(f, "  number of inner iteration = " * in_nt)
        write(f, "\n")
        write(f, "  elapse time (hh:mm:ss)    = " * et)
        write(f, "\n")
        write(f, "---------------------------- \n")
    end
end