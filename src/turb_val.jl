function set_wallpoint(nodes, bdcon, cellxmax, cellymax)
	# cal number of ponits
	nop = 0
	if Int(bdcon[1][1]) == 2 || Int(bdcon[1][1]) == 7
		nop = nop + (cellymax-1)
	end
	if Int(bdcon[2][1]) == 2 || Int(bdcon[2][1]) == 7
		nop = nop + (cellymax-1)
	end
	if Int(bdcon[3][1]) == 2 || Int(bdcon[3][1]) == 7
		nop = nop + (cellxmax-1)
	end
	if Int(bdcon[4][1]) == 2 || Int(bdcon[4][1]) == 7
		nop = nop + (cellxmax-1)
	end

	wallpoint = zeros(nop, 2) # (格子数, 座標)
	ite = 1

	# x-
	if Int(bdcon[1][1]) == 2 || Int(bdcon[1][1]) == 7
		for j in 2:cellymax
			wallpoint[ite,1] = nodes[2,j,1]
			wallpoint[ite,2] = nodes[2,j,2]
			ite = ite + 1
		end
	end
	# x+
	if Int(bdcon[2][1]) == 2 || Int(bdcon[2][1]) == 7
		for j in 2:cellymax
			wallpoint[ite,1] = nodes[2,j,1]
			wallpoint[ite,2] = nodes[2,j,2]
			ite = ite + 1
		end
	end
	# y-
	if Int(bdcon[3][1]) == 2 || Int(bdcon[3][1]) == 7
		for i in 2:cellxmax
			wallpoint[ite,1] = nodes[i,2,1]
			wallpoint[ite,2] = nodes[i,2,2]
			ite = ite + 1
		end
	end
	# y+
	if Int(bdcon[4][1]) == 2 || Int(bdcon[4][1]) == 7
		for i in 2:cellxmax
			wallpoint[ite,1] = nodes[i,2,1]
			wallpoint[ite,2] = nodes[i,2,2]
			ite = ite + 1
		end
	end
	return wallpoint, nop
end

function set_wally(wally, cellcenter, wallpoint, nop, cellxmax, cellymax)
	
	distance_list = zeros(nop,2)
	distance = 0

	pickup_point = 1

	# 最大値を1つ抽出
	for j in 2:cellymax-1
		for i in 2:cellxmax-1
			for n in 1:nop
					distance = ((cellcenter[i,j,1]-wallpoint[n,1])^2 +
								(cellcenter[i,j,2]-wallpoint[n,2])^2
								)^0.5
					distance_list[n,1] = n
					distance_list[n,2] = distance
			end

			qdistance_list = copy(distance_list)
			qdistance_list = quicksort(qdistance_list,1,nop,2,2)

			# wall_point
			wally[i,j] = qdistance_list[1,2]
			
		end
	end
	return wally
end

function cal_yplus(yplus, Qbase, wally, mu, cellxmax, cellymax)
	
	for j in 2:cellymax-1
		for i in 2:cellxmax-1
			nu = mu[i,j] / Qbase[i,j,1]
			ubar = (Qbase[i,j,2]^2 + Qbase[i,j,3]^2)^0.5
			yplus[i,j], utau = cal_wall_val_spalding(ubar, wally[i,j], nu)
			
			if isequal(yplus[i,j], NaN) == true
				println(" ")
				println(" yplus error ")
				println(Qbase[i,j,2])
				println(wally[i,j])
				println(nu)
				println(utau)
				throw(UndefVarError(:x))
			end
		end
	end

	return yplus
end

function cal_wall_val_spalding(u,y,nu)
	k = 0.4
	B = 5.5

	u_tau = 1.0
	for i in 1:20
		#=
		println("   zzz   ")
		println(u)
		println(y)
		println(nu)
		println(spalding(u_tau,u,y,nu,k,B))
		println(spalding_dash(u_tau,u,y,nu,k,B))
		println(u_tau)
		=#
		old = u_tau
		u_tau = u_tau-spalding(u_tau,u,y,nu,k,B)/spalding_dash(u_tau,u,y,nu,k,B) + 1.0e-50
		
		if abs(u_tau-old)/u_tau < 10^(-4)
			break
		end
	end

	yplus = u_tau*y/nu
	return yplus, u_tau
end

function spalding(u_tau,u,y,nu,k,B)
	t = k*u/u_tau
	F = u/u_tau-u_tau*y/nu+exp(-k*B)*(exp(-t)-1-t-t^2/2-t^3/6)
	return F
end

function spalding_dash(u_tau,u,y,nu,k,B)
	t = k*u/u_tau
	Fdash = -u/u_tau^2-y/nu+exp(-k*B)*(t/u_tau)*(-exp(-t)+1+t+t^2/2)
	return Fdash
end


# test
#=
uu = 6.64037328
yy = 0.06299976860
nunu = 1.406107133e-5
println(cal_wall_val_spalding(uu,yy,nunu))
throw(UndefVarError(:x))
=#

function wallf_Van_Driest(yplus)
	f = 1 - exp(-yplus/26)
	return f
end
