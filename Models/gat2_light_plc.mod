set V;
set O within V cross V;
set L within V cross V;

param cpt_bloc, integer >= 0;
set links_repeats within {1..cpt_bloc} cross V cross V;

param Title symbolic;
param M, integer, >=0;
param l {O};
param w {V};
param b_inf {L};
param b_sup {L};
param start symbolic;
param target symbolic;
let target := 'target_artificial';

set artificialO := setof{(pred, start) in O}(pred, target);
set newO := (O diff {(u, v) in O:v in {start}}) union artificialO;
param la{newO};

set artificialL := setof{(pred, start) in L}(pred, target);
set newL := (L diff {(u, v) in L:v in {start}}) union artificialL;
param b_inf_a {newL};
param b_sup_a {newL};

set Vi := V diff {start, target};
set reverse within V cross V;
set Flow_repeat within V cross V;
param p integer >=0;
param marge_error := 350;

param W := sum{u in V}w[u];

var x {newO}  binary;

var f {newO} integer >=0, <= W;
var i {V}   >= 0, <= 1;
var s{{start}}  >= 0, <= 1;
var t{{target}}  >= 0, <= 1;

maximize plc : sum{(u,v) in newO} x[u,v]*la[u,v] + sum{v in Vi} w[v]*i[v] + w[start];

s_start:s[start] = 1;
t_target:t[target] = 1;

#i_start:i[start] = 0;
#i_target:i[target] = 0;
no_source_reverse{(start, start_reverse) in reverse} : i[start_reverse] = 0;
no_target_reverse{(target, target_reverse) in reverse} : i[target_reverse] = 0;

#tout_le_monde_participe : sum{u in Vi}i[u] >= card(Vi)/2;
#one_edge_orientation{(u,v) in newO, (u,h) in reverse, (v,q) in reverse}: x[u,v] + x[q,h] <= 1;

visited_at_most_once{v in Vi, (v, u) in reverse} : i[v] +i[u] <= 1;
output{v in Vi} :  i[v] = sum{(v,u) in newO} x[v,u];
input{v in Vi} :  i[v] = sum{(u,v) in newO} x[u,v];

output_start :  s[start] = sum{(start,u) in newO} x[start,u];
input_start :  sum{(u,start) in newO} x[u, start] = 0;
input_target :  t[target] = sum{(u,target) in newO} x[u,target];
output_target : sum{(target, u) in newO} x[target, u] = 0;

#flow_repeat{(u,v) in Flow_repeat, (u,h) in reverse, (v,q) in reverse}: sum{(u,d) in newO} f[u,d] + sum{(h,d) in newO} f[h,d] >= sum{(k,v) in newO} f[k,v] + sum{(k,q) in newO} f[k,q];

flow_max{(u,v) in newO} : f[u,v] <= 100*W*x[u,v];
flow_output_from_source : sum{(start,u) in newO}f[start,u] = W;
subtour_elimination{v in Vi} : sum{(u,v) in newO} f[u,v] - sum{(v,u) in newO} f[v,u] = i[v]*w[v] + sum{(u,v) in newO} x[u,v]*la[u,v];