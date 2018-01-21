using DynaWAVE
using Base.Test

function test1()
    G1 = sparse([1,2,3,1],[2,1,1,3],1,3,3)
    G2 = sparse([1,2,3,1],[2,1,1,3],1,3,3)
    R = Float64[1 0 0; 0 1 0; 0 0 1]
    f = wave(G1,G2,R, 0.5)
    f == [1,2,3]
end
@test test1()

function test2()
    G1 = sparse([1,2,3,1],[2,1,1,3],1,3,3)
    G2 = sparse([1,2,3,1,4,2],[2,1,1,3,2,4],1,4,4)
    R = Float64[1 0 0 0; 0 1 0 0; 0 0 1 0]
    f = wave(G1,G2,R, 0.5)
    f == [1,2,3]
end
@test test2()

function test3()
    V = [Events([(3,5),(10,15)]), Events([(20,24),(30,40),(45,50)])]
    G1 = sparse([1,3,2,1],[2,1,1,3],vcat(V,V),3,3)
    G2 = sparse([1,3,2,1],[2,1,1,3],vcat(V,V),3,3)
    R = Float64[1 0 0; 0 1 0; 0 0 1]
    f = dynawave(G1,G2,R, 0.5)
    f == [1,2,3]
end
@test test3()
    
function test4()
    V1 = [Events([(3,5),(10,15)]), Events([(20,24),(30,40),(45,50)])]
    V2 = [Events([(3,5),(10,15)]), Events([(20,24),(30,40),(45,50)]),
          Events([(2,10),(14,20)])]
    G1 = sparse([1,3,2,1],[2,1,1,3],vcat(V1,V1),3,3)
    G2 = sparse([1,3,2,2,1,4],[2,1,4,1,3,2],vcat(V2,V2),4,4)
    R = Float64[1 0 0 0; 0 1 0 0; 0 0 1 0]
    f = dynawave(G1,G2,R, 0.5)
    f == [1,2,3]
end
@test test4()
