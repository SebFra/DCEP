set Ls := {(u, v) in L: u in {start}};
set Lt := {(u, v) in L: v in {target}};

var g {artificialL union L} binary;

maximize obj : sum{(u,v) in newL} g[u,v];
#maximize obj : sum{(u,v) in newL} g[u,v];
objective obj;

#constr_bornes_mp: p*(sum{(u,v) in L} g[u,v]) <= 300;
#distance_constr:sum{(u,v) in newO} x[u,v]*la[u,v] + sum{v in Vi} w[v]*i[v] + w[start] >= plc_inf_value;

mp_start_constraint{(u,v) in Ls} : g[u,v] <= i[v];
mp_target_constraint{(u,v) in Lt} : g[u,v] <= i[u];
mp_source_constraint{(u,v) in L diff (Ls union Lt)} : g[u,v] <= i[u];
mp_sink_constraint{(u,v) in L diff (Ls union Lt)} : g[u,v] <= i[v];
links_not_compatible{(u,v) in L, (u,q) in reverse, (v, r) in reverse}: g[u,v] + g[r, q] <= 1;

flow_distance_min_constraint_bigs{(u,v) in newL} : sum{(u,n) in newO} f[u,n] - sum{(n,v) in newO} f[n,v] + sum{(n,v) in newO} x[n,v]*la[n,v]>= (b_inf_a[u,v]- marge_error) * g[u,v] - 100*W*(1-g[u,v]);
flow_distance_max_constraint_bigs{(u,v) in newL} : sum{(u,n) in newO} f[u,n] - sum{(n,v) in newO} f[n,v] + sum{(n,v) in newO} x[n,v]*la[n,v]<= (b_sup_a[u,v] + marge_error) * g[u,v] + 100*W*(1-g[u,v]);
flow_positif_si_satisfait_constr{(u,v) in newL} : sum{(u,n) in newO} f[u,n] - sum{(n,v) in newO} f[n,v] >= (g[u,v] - 1) * W;
#test_debug1: sum{v in V} i[v] = 96;
#test_debug2: sum{(u,v) in newO} x[u,v]*la[u,v] + sum{v in Vi} w[v]*i[v] + w[start] = 175916;