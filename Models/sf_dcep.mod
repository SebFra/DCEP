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
param W := 2*sum{(u,v,c) in Edges : c = 'overlaps'}l[u,v,c];

var i{V}, >= 0, <=1;
var x{(u,v,c) in Edges_prime : c = 'overlaps'} binary;
var g{(u,v,c) in Edges : c = 'links'} binary;
var f{(u,v,c) in Edges_prime : c = 'overlaps'}, >=0, <= W;

#maximize nb_links_satisfied : sum{(u,v,'links') in Edges}g[u,v,'links'] + sum{(u,v,'overlaps') in Edges_prime} x[u,v,'overlaps']*l_st[u,v,'overlaps'];
maximize nb_links_satisfied : sum{(u,v,'links') in Edges}g[u,v,'links'];

#standart constraints
intermediate_means_just_one_predecessor{u in V} : sum{(u,v,'overlaps') in Edges_prime} x[u,v,'overlaps'] = i[u];
intermediate_input_equals_output{u in V}: sum{(u,v,'overlaps') in Edges_prime}x[u,v,'overlaps'] = sum{(v,u,'overlaps') in Edges_prime}x[v,u,'overlaps'];

source_output_constr: sum{(u,v,c) in Edges_prime : c = 'overlaps' and u = 's'} x[u,v,c] = 1;
target_output_constr: sum{(u,v,c) in Edges_prime : c = 'overlaps' and v = 't'} x[u,v,c] = 1;

visited_at_most_once{v in V, (v, u) in reverse} : i[v] +i[u] <= 1;

constraint_intermediate_node_output{u in V} : i[u] = sum{(u,v,c) in Edges_prime : c = 'overlaps'}x[u,v,c];
constraint_intermediate_node_input{u in V}  : i[u] = sum{(v,u,c) in Edges_prime : c = 'overlaps'}x[v,u,c];

no_flow_on_edge_not_on_the_path{(u,v,c) in Edges_prime : c = 'overlaps'}: f[u,v,c] <= M*x[u,v,c];

flow_initialisation:sum{(u,v,c) in Edges_prime : u = 's'}  f[u,v,c] = W;
flow_conservation{v in V}: sum{(u,v,c) in Edges_prime : c = 'overlaps'}  f[u,v,c] - sum{(v,u,c) in Edges_prime : c = 'overlaps'}  f[v,u,c] = sum{(v,u,c) in Edges_prime : c = 'overlaps'}  l_st[v,u,c]*x[v,u,c];

links_satisfied_source_must_participate{(u,v,c) in Edges : c = 'links'} : g[u,v,c] <= i[u];
links_satisfied_target_must_participate{(u,v,c) in Edges : c = 'links'} : g[u,v,c] <= i[v];

links_distance_bound_inf{(u,v,c) in Edges : c = 'links'} : sum{(u,z,'overlaps') in Edges_prime}  f[u,z,'overlaps'] - sum{(z,v,'overlaps') in Edges_prime}  f[z,v,'overlaps']  + sum{(u,z,'overlaps') in Edges_prime} l_st[u,z,'overlaps']*x[u,z,'overlaps'] >= (l[u,v,c] - marge_error)*g[u,v,c] - M*(1-g[u,v,c]);
links_distance_bound_sup{(u,v,c) in Edges : c = 'links'} : sum{(u,z,'overlaps') in Edges_prime}  f[u,z,'overlaps'] - sum{(z,v,'overlaps') in Edges_prime}  f[z,v,'overlaps']  + sum{(u,z,'overlaps') in Edges_prime} l_st[u,z,'overlaps']*x[u,z,'overlaps'] <= (l[u,v,c] + marge_error)*g[u,v,c] + M*(1-g[u,v,c]);

