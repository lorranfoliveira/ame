include("analise.jl")

using Test
using .Analise

# =======================================================================
@testset "No" begin
    ### Fields ###
    no1 = No(1, 3, 2, 5, 6, 7)
    @test no1.id == 1
    @test no1.x == 3
    @test no1.y == 2
    @test no1.fx == 5
    @test no1.fy == 6
    @test no1.mxy == 7

    no2 = No(1, 3, 2)
    @test no2.fx == 0
    @test no2.fy == 0
    @test no2.mxy == 0

    ### Erros Fields ###
    @test_throws ArgumentError No(0, 0.2, 5e9)

    ### Métodos ###
    @test Analise.distancia(No(1, -30, -10), No(2, 40, 10)) ≈ 72.80109889

    @test Analise.gl_no(No(2, 0.2, 5e9)) == [4, 5, 6]

    @test Analise.vetor_forcas(No(1, 0.2, 5e9, 3, 4, 5.6)) == [3, 4, 5.6]

end

# =======================================================================
@testset "Material" begin
    # Fields

    # id
    @test_throws ArgumentError Material(0, 0.2, 1e9)
    # E
    @test_throws ArgumentError Material(1, 0.2, -1)
    # v
    @test_throws ArgumentError Material(1, -1, 1e9)
    @test_throws ArgumentError Material(1, 1, 1e9)

    # Métodos
end

# =======================================================================
@testset "Elemento" begin
    no₁ = No(1, 3, 4)
    no₂ = No(2, 5, 6)
    mat = Material(1, 0.2, 5e9)

    ### Fields ###
    elem = Elemento(1, no₁, no₂, 10, 2, mat)   

    @test elem.id == 1
    @test elem.no₁ ≡ no₁
    @test elem.no₂ ≡ no₂
    @test elem.area == 10
    @test elem.inercia == 2
    @test elem.material ≡ mat

    ### Erros Fields ###
    # id
    @test_throws ArgumentError Elemento(0, no₁, no₂, 10, 2, Material)
    # area
    @test_throws ArgumentError Elemento(1, no₁, no₂, 0, 2, Material)
    # inercia
    @test_throws ArgumentError Elemento(1, no₁, no₂, 10, 0, Material)
end
