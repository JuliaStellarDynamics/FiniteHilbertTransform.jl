



using FiniteHilbertTransform
using OrbitalElements

const bc, M, G = 1.,1. ,1.
ψ(r::Float64)::Float64       = OrbitalElements.ψIsochrone(r,bc,M,G)
dψ(r::Float64)::Float64    = OrbitalElements.dψIsochrone(r,bc,M,G)
d2ψ(r::Float64)::Float64  = OrbitalElements.d2ψIsochrone(r,bc,M,G)
Ω₀      =    OrbitalElements.Ω₀Isochrone(bc,M,G)

n1 = -1
n2 = 2
omg = 0.0 + 0.02im
omg_nodim = omg/Ω₀  # Dimensionless frequency rescaled by Omega0
ϖ = OrbitalElements.Getϖ(omg_nodim,n1,n2,dψ,d2ψ,rmax=1000.,Ω₀=Ω₀) # Getting the rescaled frequency

println("ϖ=$ϖ")

Ku = 10

tabu,tabw,tabc,tabP = FiniteHilbertTransform.tabGLquad(Ku)
println(tabP)

FHT = FiniteHilbertTransform.LegendreFHTcreate(Ku)

println("Ku=$(FHT.Ku)")
println("tabu=$(FHT.tabu)")

# get the Legendre integration values
FiniteHilbertTransform.GettabD!(ϖ,FHT)
tabDLeg = FHT.tabDLeg
println(tabDLeg)


# [1.7935757536875392 + 1.4036549690984061im, -0.2672121231341071 - 0.7984473311775221im, -0.09183865197381542 + 0.3605951256908949im, 0.12841268709557624 - 0.1203422996242182im, -0.0847806236306678 + 0.015980864656717657im, 0.04016858917481073 + 0.016014620139537538im, -0.012898907245216886 - 0.017828308247029176im, 0.0006971307171871402 + 0.011299705804300226im, 0.0028302955525657033 - 0.005146116133935165im, -0.0026858095334769834 + 0.0014802019235564765im]

# [1.7935757537288612 + 1.403654969105378im, -0.267212123157653 - 0.7984473311962709im, -0.09183865196684303 + 0.36059512570936236im, 0.12841268709788553 - 0.12034229963638196im, -0.08478062363577567 + 0.015980864662294734im, 0.04016858917918258 + 0.01601462013831581im, -0.012898907247798562 - 0.017828308247788426im, 0.0006971307182353672 + 0.011299705805483163im, 0.0028302955524146466 - 0.005146116134844375im, -0.0026858095336833037 + 0.0014802019240494299im]
