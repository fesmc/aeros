using NCDatasets, Statistics, Printf
ds = NCDataset(joinpath(@__DIR__,"output","speedy_omega_T21L8.nc"))
lat = ds["lat"][:]; sig = ds["sigma"][:]
omega = ds["omega"][:,:]; rh = ds["rh"][:,:]; heat = ds["heat"][:,:]
pc = ds["precip_convective"][:]; pl = ds["precip_largescale"][:]
nlat, nlev = size(omega)

# hemispheric symmetrization onto |lat| bins (lat is ascending -85..85)
absl = abs.(lat)
ubins = sort(unique(round.(absl; digits=1)))
sym(field2d) = [mean(field2d[findall(x->isapprox(round(abs(x);digits=1),b), lat), k]) for b in ubins, k in 1:nlev]
symv(v) = [mean(v[findall(x->isapprox(round(abs(x);digits=1),b), lat)]) for b in ubins]
O = sym(omega); R = sym(rh); H = sym(heat); PC = symv(pc); PL = symv(pl)

println("=== hemispherically-symmetrized |lat| profiles ===")
println("sigma levels: ", round.(sig;digits=3))
@printf("%6s | %s | precip(mm/d)\n","|lat|","omega(hPa/day) by sigma")
for (i,b) in enumerate(ubins)
    @printf("%6.1f |", b)
    for k in 1:nlev; @printf(" %6.2f", O[i,k]); end
    @printf(" | %5.2f\n", PC[i]+PL[i])
end
println("\n=== RH(|lat|,sigma) % ===")
for (i,b) in enumerate(ubins)
    @printf("%6.1f |", b); for k in 1:nlev; @printf(" %5.1f", R[i,k]); end; println()
end

# headline numbers (symmetrized)
trop = findall(b->b<=25, ubins)
# tropical peak ascent (most negative over levels 2:nlev-1)
casc = [minimum(O[i,2:nlev-1]) for i in 1:length(ubins)]
it = trop[argmin(casc[trop])]
@printf("\nTROP peak ascent (sym): %.2f hPa/day at |lat|=%.1f\n", casc[it], ubins[it])
# subtropical descent 15-35
cdesc = [maximum(O[i,2:nlev-1]) for i in 1:length(ubins)]
sb = findall(b->15<=b<=35, ubins)
id = sb[argmax(cdesc[sb])]
@printf("SUBTROP descent (sym): %.2f hPa/day at |lat|=%.1f\n", cdesc[id], ubins[id])
# equatorial omega (near 0-3 deg) mid-trop
eq = argmin(ubins)
@printf("EQUATOR omega mid-trop (sig~0.44): %.2f hPa/day\n", O[eq,4])
# subtropical free-trop RH minimum band 25-40, mid trop (sig 0.31-0.56)
sb2 = findall(b->25<=b<=40, ubins)
@printf("SUBTROP free-trop RH min (25-40, sig0.31-0.56): %.1f %%\n", minimum(R[sb2,3:5]))
@printf("SUBTROP free-trop RH mean (15-35, sig~0.44): %.1f %%\n", mean(R[findall(b->15<=b<=35,ubins),4]))
# deep-tropics precip: equatorial vs off-equatorial
@printf("PRECIP: equator(|lat|~2.8)=%.2f  off-eq peak=%.2f at |lat|=%.1f\n",
        PC[eq]+PL[eq], maximum(PC.+PL), ubins[argmax(PC.+PL)])
close(ds)
