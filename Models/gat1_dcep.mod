param Title symbolic;
param start symbolic;
set V;
set reverse within V cross V;
set Edges within V cross V cross {'overlaps', 'links'};

param w{V} integer >= 0;
param l{Edges} integer;
param marge_error := 350;
param M := 100*sum{(u,v,c) in Edges}l[u,v,c];
param W := 10*sum{(u,v,c) in Edges : c = 'overlaps'}l[u,v,c];

var i{V}, >= 0, <=1;
var s{V}, >= 0, <=1;
var t{V}, >= 0, <=1;
var x{(u,v,c) in Edges : c = 'overlaps'} binary;
var g{(u,v,c) in Edges : c = 'links'} binary;
var f{(u,v,c) in Edges : c = 'overlaps'}, >=0, <= W;

#maximize nb_links_satisfied : M*sum{(u,v,'links') in Edges}g[u,v,'links'] + sum{(u,v,'overlaps') in Edges}f[u,v,'overlaps'];
maximize nb_links_satisfied : sum{(u,v,'links') in Edges}g[u,v,'links'];

#standart constraints
intermediate_means_just_one_predecessor{u in V} : sum{(u,v,'overlaps') in Edges} x[u,v,'overlaps'] <= 1;

only_one_source : sum{v in V}s[v] = 1;
only_one_target : sum{v in V}t[v] = 1;
output{v in V} :  i[v] + s[v] = sum{(v,u,'overlaps') in Edges} x[v,u, 'overlaps'];
input{v in V} :  t[v] + i[v] = sum{(u,v,'overlaps') in Edges} x[u,v, 'overlaps'];

visited_at_most_once{v in V, (v, u) in reverse} : i[v] + s[v] + t[v]+i[u] + s[u] + t[u] <= 1;

no_flow_on_edge_not_on_the_path{(u,v,c) in Edges : c = 'overlaps'}: f[u,v,c] <= W*x[u,v,c];

#flow_conservation{v in V}: sum{(u,v,c) in Edges : c = 'overlaps'}  f[u,v,c] - sum{(v,u,c) in Edges : c = 'overlaps'}  f[v,u,c] = sum{(v,u,c) in Edges : c = 'overlaps'}  l[v,u,c]*x[v,u,c];
flow_output_from_source{v in V} : W*s[v] <= sum{(v,u,'overlaps') in Edges}f[v,u,'overlaps'] + sum{(v,u, 'overlaps') in Edges} x[v,u,'overlaps']*l[v,u,'overlaps'];

subtour_elimination{v in V} : sum{(u,v,'overlaps') in Edges} f[u,v,'overlaps'] - sum{(v,u,'overlaps') in Edges} f[v,u,'overlaps']  >=  sum{(v,u, 'overlaps') in Edges} x[v,u,'overlaps']*l[v,u,'overlaps'] - W*(s[v]);

links_satisfied_source_must_participate{(u,v,c) in Edges : c = 'links'} : g[u,v,c] <= s[u]+i[u]+t[u];
links_satisfied_target_must_participate{(u,v,c) in Edges : c = 'links'} : g[u,v,c] <= s[v]+i[v]+t[v];

links_distance_bound_inf{(u,v,c) in Edges : c = 'links'} : sum{(u,z,'overlaps') in Edges}  f[u,z,'overlaps'] - sum{(z,v,'overlaps') in Edges}  f[z,v,'overlaps']  + sum{(u,z,'overlaps') in Edges} l[u,z,'overlaps']*x[u,z,'overlaps'] >= (l[u,v,c] - marge_error)*g[u,v,c] - M*(1-g[u,v,c]);
links_distance_bound_sup{(u,v,c) in Edges : c = 'links'} : sum{(u,z,'overlaps') in Edges}  f[u,z,'overlaps'] - sum{(z,v,'overlaps') in Edges}  f[z,v,'overlaps']  + sum{(u,z,'overlaps') in Edges} l[u,z,'overlaps']*x[u,z,'overlaps'] <= (l[u,v,c] + marge_error)*g[u,v,c] + M*(1-g[u,v,c]);
