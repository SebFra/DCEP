param Title symbolic;
param start symbolic;
set V;
set reverse within V cross V;
set Edges within V cross V cross {'overlaps', 'links'};
set V_prime := V union {'s', 't'};
set Edges_prime dimen 3;

param w{V} integer >= 0;
param l{Edges} integer;
param l_st{Edges_prime} integer;
param marge_error := 350;
param M := 100*sum{(u,v,c) in Edges}l[u,v,c];
param W := sum{(u,v,c) in Edges : c = 'overlaps'}l[u,v,c];

var x{(u,v,c) in Edges_prime : c = 'overlaps'} binary;
var i{V}, binary;
var y{V_prime} integer >= 0, <=W;
var g{(u,v,c) in Edges : c = 'links'} binary;

maximize nb_links_satisfied : M*sum{(u,v,'links') in Edges}g[u,v,'links'] + sum{v in V}y[v];

#standart constraints
intermediate_means_just_one_predecessor{u in V} : sum{(u,v,'overlaps') in Edges_prime} x[u,v,'overlaps'] = i[u];
intermediate_input_equals_output{u in V}: sum{(u,v,'overlaps') in Edges_prime}x[u,v,'overlaps'] = sum{(v,u,'overlaps') in Edges_prime}x[v,u,'overlaps'];

initialisation_y_s : y['s'] = W;
source_output_constr: sum{(u,v,c) in Edges_prime : u = 's'} x[u,v,c] = 1;
target_output_constr: sum{(u,v,c) in Edges_prime : v = 't'} x[u,v,c] = 1;

visited_at_most_once{v in V, (v, u) in reverse} : i[v] +i[u] <= 1;

subtour_elimination{(u,v,'overlaps') in Edges_prime}: y[v] <= y[u] - x[u,v,'overlaps']*l_st[u,v,'overlaps'] + (1-x[u,v,'overlaps'])*M;
label_cancelation{v in V}: y[v] <= M*sum{(u,v,'overlaps') in Edges_prime} x[u,v,'overlaps'];

links_satisfied_source_must_participate{(u,v,c) in Edges : c = 'links'} : g[u,v,c] <= y[u];
links_satisfied_target_must_participate{(u,v,c) in Edges : c = 'links'} : g[u,v,c] <= y[v];

links_distance_bound_inf{(u,v,c) in Edges : c = 'links'} : y[u] - y[v] >= (l[u,v,c] - marge_error)*g[u,v,c] - M*(1-g[u,v,c]);
links_distance_bound_sup{(u,v,c) in Edges : c = 'links'} : y[u] - y[v]  <= (l[u,v,c] + marge_error)*g[u,v,c] + M*(1-g[u,v,c]);

