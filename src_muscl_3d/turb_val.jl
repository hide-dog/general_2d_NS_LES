function set_wally(nodes, bdcon, wally, cellcenter, cellxmax, cellymax, cellzmax, icell)
	# serch wall
	swith_wall = zeros(6)   # x-, x+, y-, y+, z-, z+を検索する
							# 1:壁
	for i in 1:6
		if Int(bdcon[i][1]) == 2 || Int(bdcon[i][1]) == 7
			swith_wall[i] = 1
		end
	end

	# cal number of ponits
	nop = 0
	if swith_wall[1] == 1
		nop = nop + (cellymax-1)*(cellzmax-1)
	end
	if swith_wall[2] == 1
		nop = nop + (cellymax-1)*(cellzmax-1)
	end
	if swith_wall[3] == 1
		nop = nop + (cellzmax-1)*(cellxmax-1)
	end
	if swith_wall[4] == 1
		nop = nop + (cellzmax-1)*(cellxmax-1)
	end
	if swith_wall[5] == 1
		nop = nop + (cellxmax-1)*(cellzmax-1)
	end
	if swith_wall[6] == 1
		nop = nop + (cellxmax-1)*(cellzmax-1)
	end

	wallpoint = zeros(nop, 3) # (格子数, 座標)
	ite = 1

	# x-
	if swith_wall[1] == 1
		for k in 1+icell:cellzmax-icell
			for j in 1+icell:cellymax-icell
				wallpoint[ite,1] = nodes[icell+1,j,k,1]
				wallpoint[ite,2] = nodes[icell+1,j,k,2]
				wallpoint[ite,3] = nodes[icell+1,j,k,3]
				ite = ite + 1
			end
		end
	end
	# x+
	if swith_wall[2] == 1
		for k in 1+icell:cellzmax-icell
			for j in 1+icell:cellymax-icell
				wallpoint[ite,1] = nodes[cellxmax-icell+1,j,k,1]
				wallpoint[ite,2] = nodes[cellxmax-icell+1,j,k,2]
				wallpoint[ite,3] = nodes[cellxmax-icell+1,j,k,3]
				ite = ite + 1
			end
		end
	end
	# y-
	if swith_wall[3] == 1
		for k in 1+icell:cellzmax-icell
			for i in 1+icell:cellxmax-icell
				wallpoint[ite,1] = nodes[i,icell+1,k,1]
				wallpoint[ite,2] = nodes[i,icell+1,k,2]
				wallpoint[ite,3] = nodes[i,icell+1,k,3]
				ite = ite + 1
			end
		end
	end
	# y+
	if swith_wall[4] == 1
		for k in 1+icell:cellzmax-icell
			for i in 1+icell:cellxmax-icell
				wallpoint[ite,1] = nodes[i,cellymax-icell+1,k,1]
				wallpoint[ite,2] = nodes[i,cellymax-icell+1,k,2]
				wallpoint[ite,3] = nodes[i,cellymax-icell+1,k,3]
				ite = ite + 1
			end
		end
	end
	# z-
	if swith_wall[5] == 1
		for j in 1+icell:cellymax-icell
			for i in 1+icell:cellxmax-icell
				wallpoint[ite,1] = nodes[i,j,icell+1,1]
				wallpoint[ite,2] = nodes[i,j,icell+1,2]
				wallpoint[ite,3] = nodes[i,j,icell+1,3]
				ite = ite + 1
			end
		end
	end
	# z+
	if swith_wall[6] == 1
		for j in 1+icell:cellymax-icell
			for i in 1+icell:cellxmax-icell
				wallpoint[ite,1] = nodes[i,j,cellzmax-icell+1,1]
				wallpoint[ite,2] = nodes[i,j,cellzmax-icell+1,2]
				wallpoint[ite,3] = nodes[i,j,cellzmax-icell+1,3]
				ite = ite + 1
			end
		end
	end



	# 距離格納用
	distance_list = zeros(nop, 3)
	distance = 0

	# 最大値を2つ抽出
	for k in 1+icell:cellymax-icell
		for j in 1+icell:cellymax-icell
			for i in 1+icell:cellxmax-icell
				for n in 1:nop
					distance = ((cellcenter[i,j,k,1]-wallpoint[n,1])^2 +
								(cellcenter[i,j,k,2]-wallpoint[n,2])^2 +
								(cellcenter[i,j,k,3]-wallpoint[n,3])^2
								)^0.5
					distance_list[n,1] = n
					distance_list[n,2] = distance
				end

				qdistance_list = copy(distance_list)
				qdistance_list = quicksort(qdistance_list,1,nop,2,2)

				# wall_point
				point1 = Int(qdistance_list[1,1])
				point2 = Int(qdistance_list[2,1])
				
				dis_bottom = ((wallpoint[point1, 1] - wallpoint[point2, 1])^2 +
							  (wallpoint[point1, 2] - wallpoint[point2, 2])^2 +
							  (wallpoint[point1, 3] - wallpoint[point2, 3])^2)^0.5
				
				# O型格子の関係で同じ点が二つあるため
				if dis_bottom == 0.0
					point1 = Int(qdistance_list[1,1])
					point2 = Int(qdistance_list[3,1])
				
					dis_bottom = ((wallpoint[point1, 1] - wallpoint[point2, 1])^2 +
				   			      (wallpoint[point1, 2] - wallpoint[point2, 2])^2 +
							      (wallpoint[point1, 3] - wallpoint[point2, 3])^2)^0.5
				end

				dis1 = distance_list[point1, 2]
				dis2 = distance_list[point2, 2]
				
				# ヘロンの公式から距離を算出
				# 正確な平面への投影はしてないので、注意
				s = (dis_bottom + dis1 + dis2)*0.5
				surface = (s*(s-dis_bottom)*(s-dis1)*(s-dis2))^0.5
				wally[i,j,k] = surface/dis_bottom * 2
			end
		end
	end

	# 仮想セルへの代入
	for k in 1+icell:cellzmax-icell
		for j in 1:icell
			for i in 1+icell:cellxmax-icell
				wally[i,j,k] = wally[i,icell+1,k]
				wally[i,cellymax-j+1,k] = wally[i,cellymax-icell,k]
			end
		end
	end
	for k in 1+icell:cellzmax-icell
		for j in 1+icell:cellymax-icell
			for i in 1:icell
				wally[i,j,k] = wally[icell+1,j,k]
				wally[cellxmax-i+1,j,k] = wally[cellxmax-icell,j,k]
			end
		end
	end
	for k in 1:icell
		for j in 1+icell:cellymax-icell
			for i in 1+icell:cellxmax-icell
				wally[i,j,k] = wally[i,j,icell+1]
				wally[i,j,cellzmax-k+1] = wally[i,j,cellzmax-icell]
			end
		end
	end

	return wally, swith_wall
end

function cal_yplus(yplus, Qbase, wally, swith_wall, mu, cellxmax, cellymax, cellzmax, vecAx, vecAy, vecAz, volume, icell)
	tvA = zeros(3)
	
	# 壁に並行な方向の速度を求める
	# 速度ベクトル内壁に垂直な成分は、ux = vec{u}*vecAx
	# u_purpose = (u^2+v^2+w^2 - ux^2)^2
	for k in 1+icell:cellzmax-icell
		for j in 1+icell:cellymax-icell
			for i in 1+icell:cellzmax-icell
				if swith_wall[1] == 1 || swith_wall[2] == 1
					tvA[1] = 0.5*(vecAx[i-1,j,k,1] + vecAx[i,j,k,1])
					tvA[2] = 0.5*(vecAx[i-1,j,k,2] + vecAx[i,j,k,2])
					tvA[3] = 0.5*(vecAx[i-1,j,k,3] + vecAx[i,j,k,3])
				elseif swith_wall[3] == 1 || swith_wall[4] == 1
					tvA[1] = 0.5*(vecAy[i,j-1,k,1] + vecAy[i,j,k,1])
					tvA[2] = 0.5*(vecAy[i,j-1,k,2] + vecAy[i,j,k,2])
					tvA[3] = 0.5*(vecAy[i,j-1,k,3] + vecAy[i,j,k,3])
				elseif swith_wall[5] == 1 || swith_wall[6] == 1
					tvA[1] = 0.5*(vecAz[i,j,k-1,1] + vecAz[i,j,k,1])
					tvA[2] = 0.5*(vecAz[i,j,k-1,2] + vecAz[i,j,k,2])
					tvA[3] = 0.5*(vecAz[i,j,k-1,3] + vecAz[i,j,k,3])
				end

				sqA = (tvA[1]^2 + tvA[2]^2 + tvA[3]^2)^0.5
				tvA[1] = tvA[1] / sqA
				tvA[2] = tvA[2] / sqA
				tvA[3] = tvA[3] / sqA

				nu = mu[i,j,k] / Qbase[i,j,k,1]
				u = (tvA[1]*Qbase[i,j,k,2] + tvA[2]*Qbase[i,j,k,3] + tvA[3]*Qbase[i,j,k,4])^0.5
				ubar = (Qbase[i,j,k,2]^2 + Qbase[i,j,k,3]^2 + Qbase[i,j,k,4]^2 - u^2)^0.5
				yplus[i,j,k], up = cal_wall_val_spalding(ubar, wally[i,j,k], nu)
			end
		end
	end

	return yplus
end

function cal_wall_val_spalding(u,y,nu)
	# spalding 定数
	k = 0.4
	B = 5.5

	# ニュートン法
	up = 200.0
	for i in 1:100
		old = up
		
		up = up - spalding(up,u,y,nu,k,B)/spalding_dash(up,u,y,nu,k,B)

		if abs(up-old)/up < 10^(-4)
			break
		end
	end

	yplus = u*y/nu/up
	return yplus, up
end

function spalding(up,u,y,nu,k,B)
	t = k*up
	F = - 1/up* u*y/nu + up + exp(-k*B)*(exp(t)-1-t-t^2/2-t^3/6)
	return F
end

function spalding_dash(up,u,y,nu,k,B)
	t = k*up
	Fdash = 1/up^2* u*y/nu + 1 + exp(-k*B)*(k*exp(t)-k-k^2*up-k^3*up^2/2)
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
