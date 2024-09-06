
function tbf90(I01p, T, incAng, k, a, g, n, eps, dz)
    cos1 = cos(incAng / 180.0 * π)
    Jup = zeros(eltype(I01p), n + 1)
    Jd = zeros(eltype(I01p), n + 1)

    intDo = zeros(eltype(I01p), n + 1)
    intUp = zeros(eltype(I01p), n + 1)
    
    # these temperatures should be defined at the boundaries between cells
    # there are n I0 and n+1 I1 values
    for k1 in 1:n+1
        Jup[k1] = (1 - a[k1]) * T[k1] + a[k1] * (0.5 * (I01p[2*k1-1] + I01p[2*k1+1]) * g[k1] * cos1 + I01p[2*k1])
        Jd[k1] = (1 - a[k1]) * T[k1] + a[k1] * (-0.5 * (I01p[2*k1-1] + I01p[2*k1+1]) * g[k1] * cos1 + I01p[2*k1])
    end

    intDo[1] = k[1] * dz / 2
    intUp[n+1] = k[n+1] * dz / 2

    for k1 in 1:n
        intUp[n+1-k1] = intUp[n-k1+2] + 0.5 * (k[n+1-k1] + k[n+1-k1+1]) * dz
    end

    for k1 in 2:n+1
        intDo[k1] = intDo[k1-1] + 0.5 * (k[k1] + k[k1-1]) * dz
    end
    sumJD = zero(eltype(I01p))
    sumJU = zero(eltype(I01p))
    
    for i in 1:n+1
        sumJD += Jd[i] * exp(-intDo[i] / cos1) * k[i] * dz / cos1
        sumJU += Jup[i] * exp(-intUp[i] / cos1) * k[i] * dz / cos1
    end

    Tb = (1 - eps) * exp(-intUp[1] / cos1) * sumJD + sumJU + eps * exp(-intUp[1] / cos1) * T[1]
    return Tb
end

function SetEddington1D(T, k, a, g, n, eps, dz, Ts)
    Abig = zeros(eltype(T), 2*n+3, 2*n+3)
    B = zeros(eltype(T), 2*n+3)
    # dI0/dz=-k*(1-a*g)*I1      deq  
    # dI1/dz=-3*k*(1-a)*(I0-T)  deq 
    # I0,0+(2-eps)/(3*eps)*(I1,0+I1,1)=Ts bc
    # the first component of the solution is I1,1, then I0,1 and so on
    # J=(1-a)*T+a*(I0+g*I1*cos(ang))  J
    # Io,n -1/3 [I1,n+I1,n+1]= 2.7 top boundary
    # In this discretization, temperatures and scattering properties should be defined at the centers of cells
    for i in 0:n
        Abig[2*i+2, 2*i+2] = 3 * k[i+1] * (1 - a[i+1]) * dz  # equation for I1
        Abig[2*i+2, 2*i+3] = 1
        Abig[2*i+2, 2*i+1] = -1
        B[2*i+2] = 3 * k[i+1] * (1 - a[i+1]) * T[i+1] * dz

        if 2*i+3 < 2*n+3
            km = 0.5 * (k[i+1] + k[i+2])
            am = 0.5 * (a[i+1] + a[i+2])
            gm = 0.5 * (g[i+1] + g[i+2])
            Abig[2*i+3, 2*i+3] = km * (1 - am * gm) * dz # equation for I0
            Abig[2*i+3, 2*i+4] = 1
            Abig[2*i+3, 2*i+2] = -1
            B[2*i+3] = 0.0
        end
    end

    Abig[2*n+3, 2*n+2] = 1.0
    Abig[2*n+3, 2*n+3] = -1.0 / 3.0
    Abig[2*n+3, 2*n+1] = -1.0 / 3.0
    B[2*n+3] = 2.7

    Abig[1, 1] = (2 - eps) / (3 * eps)
    Abig[1, 2] = 1.0
    Abig[1, 3] = (2 - eps) / (3 * eps)
    B[1] = Ts
    return Abig, B 
end

function eddingtonTb(temp1, kext1, kscat1, asym1, n, eps, dz, incAng)
    temp_mid=temp1;
    A,B=SetEddington1D(temp_mid, kext1, kscat1, asym1, n, eps, dz, temp1[1]);
    I01p=A\B;
    Tb_jl=tbf90(I01p, temp1, incAng, kext1, kscat1, asym1, n, eps, dz);
    return Tb_jl
end