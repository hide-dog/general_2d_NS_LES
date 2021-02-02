function set_wally(nodes, bdcon, wally, cellcenter, cellxmax, cellymax)
	# serch wall
	swith_wall = zeros(4)   # x-, x+, y-, y+を検索する
							# 1:壁
	for i in 1:4
		if Int(bdcon[i][1]) == 2 || Int(bdcon[i][1]) == 7
			swith_wall[i] = 1
		end
	end

	# cal number of ponits
	nop = 0
	if swith_wall[1] == 1
		nop = nop + (cellymax-1)
	end
	if swith_wall[2] == 1
		nop = nop + (cellymax-1)
	end
	if swith_wall[3] == 1
		nop = nop + (cellxmax-1)
	end
	if swith_wall[4] == 1
		nop = nop + (cellxmax-1)
	end

	wallpoint = zeros(nop, 2) # (格子数, 座標)
	ite = 1

	# x-
	if swith_wall[1] == 1
		for j in 2:cellymax-1
			wallpoint[ite,1] = nodes[2,j,1]
			wallpoint[ite,2] = nodes[2,j,2]
			ite = ite + 1
		end
	end
	# x+
	if swith_wall[2] == 1
		for j in 2:cellymax-1
			wallpoint[ite,1] = nodes[2,j,1]
			wallpoint[ite,2] = nodes[2,j,2]
			ite = ite + 1
		end
	end
	# y-
	if swith_wall[3] == 1
		for i in 2:cellxmax-1
			wallpoint[ite,1] = nodes[i,2,1]
			wallpoint[ite,2] = nodes[i,2,2]
			ite = ite + 1
		end
	end
	# y+
	if swith_wall[4] == 1
		for i in 2:cellxmax-1
			wallpoint[ite,1] = nodes[i,2,1]
			wallpoint[ite,2] = nodes[i,2,2]
			ite = ite + 1
		end
	end

	# 距離格納用
	distance_list = zeros(nop,2)
	distance = 0

	# 最大値を2つ抽出
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
			point1 = Int(qdistance_list[1,1])
			point2 = Int(qdistance_list[2,1])
			
			dis_bottom = ((wallpoint[point1, 1] - wallpoint[point2, 1])^2 + (wallpoint[point1, 2] - wallpoint[point2, 2])^2)^0.5
			
			# O型格子の関係で同じ点が二つあるため
			if dis_bottom == 0.0
				point1 = Int(qdistance_list[1,1])
				point2 = Int(qdistance_list[3,1])
			
				dis_bottom = ((wallpoint[point1, 1] - wallpoint[point2, 1])^2 + (wallpoint[point1, 2] - wallpoint[point2, 2])^2)^0.5
			end

			dis1 = distance_list[point1, 2]
			dis2 = distance_list[point2, 2]
			
			# ヘロンの公式から距離を算出
			s = (dis_bottom + dis1 + dis2)*0.5
			surface = (s*(s-dis_bottom)*(s-dis1)*(s-dis2))^0.5
			wally[i,j] = surface/dis_bottom * 2
		end
	end
	return wally, swith_wall
end

function cal_yplus(yplus, Qbase, wally, swith_wall, mu, cellxmax, cellymax, vecAx, vecAy)
	
	if swith_wall[1] == 1
		for j in 2:cellymax-1
			ite = [2, 3]
			for i in ite
				Axx = 0.5*(vecAx[i-1,j,1] + vecAx[i,j,1])
				Axy = 0.5*(vecAx[i-1,j,2] + vecAx[i,j,2])
				nu = mu[i,j] / Qbase[i,j,1]
				ubar = abs(Axx*Qbase[i,j,2] + Axy*Qbase[i,j,3])
				yplus[i,j], utau = cal_wall_val_spalding(ubar, wally[i,j], nu)
			end
		end
	end
	if swith_wall[2] == 1
		for j in 2:cellymax-1
			ite = [cellxmax-2, cellxmax-1]
			for i in ite
				Axx = 0.5*(vecAx[i-1,j,1] + vecAx[i,j,1])
				Axy = 0.5*(vecAx[i-1,j,2] + vecAx[i,j,2])
				nu = mu[i,j] / Qbase[i,j,1]
				ubar = abs(Axx*Qbase[i,j,2] + Axy*Qbase[i,j,3])
				yplus[i,j], utau = cal_wall_val_spalding(ubar, wally[i,j], nu)
			end
		end
	end
	if swith_wall[3] == 1
		for i in 2:cellxmax-1
			ite = [2, 3]
			for j in ite
				Ayx = 0.5*(vecAy[i,j-1,1] + vecAx[i,j,1])
				Ayy = 0.5*(vecAy[i,j-1,2] + vecAx[i,j,2])
				nu = mu[i,j] / Qbase[i,j,1]
				ubar = abs(Ayx*Qbase[i,j,2] + Ayy*Qbase[i,j,3])
				#yplus[i,j], utau = cal_wall_val_spalding(ubar, wally[i,j], nu)
				#yplus[i,j] = 30.0
			end
		end
	end
	if swith_wall[4] == 1
		for i in 2:cellxmax-1
			ite = [cellymax-2, cellymax-1]
			for j in ite
				Ayx = 0.5*(vecAy[i,j-1,1] + vecAx[i,j,1])
				Ayy = 0.5*(vecAy[i,j-1,2] + vecAx[i,j,2])
				nu = mu[i,j] / Qbase[i,j,1]
				ubar = abs(Ayx*Qbase[i,j,2] + Ayy*Qbase[i,j,3])
				yplus[i,j], utau = cal_wall_val_spalding(ubar, wally[i,j], nu)
			end
		end
	end

	return yplus
end


function cal_wall_val_spalding(u,y,nu)
	# spalding 定数
	k = 0.4
	B = 5.5

	# 二分法
	u_taup = 500
	u_taum = -10
	u_tau = 1.0
	for i in 1:20
		old = u_tau
		
		# 減速項nuの計算ループ
		for j in 1:5

			temp_u_tau = u_tau - nu * spalding(u_tau,u,y,nu,k,B)/spalding_dash(u_tau,u,y,nu,k,B)

			if abs(spalding(temp_u_tau,u,y,nu,k,B)) < (1-nu*sigma) * abs(spalding(temp_u_tau,u,y,nu,k,B))
				u_tau = temp_u_tau
				break
			else
				nu = nu/2
			end
		end

		if abs(u_tau-old)/u_tau < 10^(-4)
			break
		end
	end

	yplus = u_tau*y/nu
	return yplus, u_tau
end


function cal_wall_val_spalding_newton(u,y,nu)
	# spalding 定数
	k = 0.4
	B = 5.5

	# 減速ニュートン法の定数
	nu = 1
	sigma = 0.25 # <1

	# 減速ニュートン法
	u_tau = 1
	for i in 1:20
		old = u_tau
		
		# 減速項nuの計算ループ
		for j in 1:5

			temp_u_tau = u_tau - nu * spalding(u_tau,u,y,nu,k,B)/spalding_dash(u_tau,u,y,nu,k,B)

			if abs(spalding(temp_u_tau,u,y,nu,k,B)) < (1-nu*sigma) * abs(spalding(temp_u_tau,u,y,nu,k,B))
				u_tau = temp_u_tau
				break
			else
				nu = nu/2
			end
		end

		if abs(u_tau-old)/u_tau < 10^(-4)
			break
		end
	end

	yplus = u_tau*y/nu
	return yplus, u_tau
end

function spalding(u_tau,u,y,nu,k,B)
	t = k*u/u_tau
	F = u/u_tau-u_tau*y/nu + exp(-k*B)*(exp(t)-1-t-t^2/2-t^3/6)
	return F
end

function spalding_dash(u_tau,u,y,nu,k,B)
	t = k*u/u_tau
	Fdash = -u/u_tau^2-y/nu + exp(-k*B)*(t/u_tau)*(-exp(t)+1+t+t^2/2)
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
