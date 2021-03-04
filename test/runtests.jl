using SPHKernels, Test

@testset "SPH Kernels" begin

    @testset "Cubic Spline" begin

        @testset "1D value" begin
            k = Cubic()
            # < 0.5
            d = kernel_value_1D(k, 0.4, 0.5)
            @test d ≈ 0.2826666666666666       
            # < 1.0
            d = kernel_value_1D(k, 0.5, 0.5)
            @test d ≈ 0.16666666666666666
            # > 1.0
            d = kernel_value_1D(k, 1.5, 0.5)
            @test d == 0.0
        end

        @testset "1D derivative" begin
            k = Cubic()
            # < 0.5
            d = kernel_deriv_1D(k, 0.4, 0.5)
            @test d ≈ -0.6399999999999999    
            # < 1.0
            d = kernel_deriv_1D(k, 0.5, 0.5)
            @test d ≈ -0.5
            # > 1.0
            d = kernel_deriv_1D(k, 1.5, 0.5)
            @test d == 0.0
        end

        @testset "1D bias correction" begin
            k = Cubic()
            d = bias_correction_1D(k, 1.0, 0.0, 0.0)
            @test d ≈  1.0  
        end

        @testset "2D value" begin
            k = Cubic()
            # < 0.5
            d = kernel_value_2D(k, 0.4, 0.5)
            @test d ≈ 0.19280484534561032       
            # < 1.0
            d = kernel_value_2D(k, 0.5, 0.5)
            @test d ≈ 0.11368210220849667
            # > 1.0
            d = kernel_value_2D(k, 1.5, 0.5)
            @test d == 0.0
        end

        @testset "2D derivative" begin
            k = Cubic()
            # < 0.5
            d = kernel_deriv_2D(k, 0.4, 0.5)
            @test d ≈ -0.4365392724806272     
            # < 1.0
            d = kernel_deriv_2D(k, 0.5, 0.5)
            @test d ≈ -0.34104630662549
            # > 1.0
            d = kernel_deriv_2D(k, 1.5, 0.5)
            @test d == 0.0
        end

        @testset "2D bias correction" begin
            k = Cubic()
            d = bias_correction_2D(k, 1.0, 0.0, 0.0)
            @test d ≈  1.0  
        end
        
        @testset "3D value" begin
            k = Cubic()
            # < 0.5
            d = kernel_value_3D(k, 0.4, 0.5)
            @test d ≈ 0.13496339174192723
            # < 1.0
            d = kernel_value_3D(k, 0.5, 0.5)
            @test d ≈ 0.07957747154594767
            # > 1.0
            d = kernel_value_3D(k, 1.5, 0.5)
            @test d == 0.0
        end

        @testset "3D derivative" begin
            k = Cubic()
            # < 0.5
            d = kernel_deriv_3D(k, 0.4, 0.5)
            @test d ≈ -0.30557749073643903    
            # < 1.0
            d = kernel_deriv_3D(k, 0.5, 0.5)
            @test d ≈ -0.238732414637843
            # > 1.0
            d = kernel_deriv_3D(k, 1.5, 0.5)
            @test d == 0.0
        end

        @testset "3D bias correction" begin
            k = Cubic()
            d = bias_correction_3D(k, 1.0, 0.0, 0.0)
            @test d ≈  1.0  
        end
    end
    
    @testset "Quintic Spline" begin
        
        @testset "1D value" begin
            k = Quintic()
            # < 1/3
            d = kernel_value_1D(k, 0.3, 0.5)
            @test d ≈ 0.38972624999999994
            # < 2/3
            d = kernel_value_1D(k, 0.5, 0.5)
            @test d ≈ 0.092578125
            # < 1.0
            d = kernel_value_1D(k, 0.8, 0.5)
            @test d ≈ 0.0009719999999999988
            # > 1.0
            d = kernel_value_1D(k, 1.5, 0.5)
            @test d == 0.0
        end

        @testset "1D derivative" begin
            k = Quintic()
            # < 1/3
            d = kernel_deriv_1D(k, 0.3, 0.5)
            @test d ≈ -0.99984375
            # < 2/3
            d = kernel_deriv_1D(k, 0.5, 0.5)
            @test d ≈ -114.345703125
            # < 1.0
            d = kernel_deriv_1D(k, 0.8, 0.5)
            @test d ≈ -0.012149999999999987
            # > 1.0
            d = kernel_deriv_1D(k, 1.5, 0.5)
            @test d == 0.0
        end

        @testset "1D bias correction" begin
            k = Quintic()
            d = bias_correction_1D(k, 1.0, 0.0, 0.0)
            @test d ≈  1.0  
        end

        @testset "2D value" begin
            k = Quintic()
            # < 1/3
            d = kernel_value_2D(k, 0.3, 0.5)
            @test d ≈ 0.32700352517410614
            # < 2/3
            d = kernel_value_2D(k, 0.5, 0.5)
            @test d ≈ 0.07767855829318415
            # < 1.0
            d = kernel_value_2D(k, 0.8, 0.5)
            @test d ≈ 0.0008155658657050455
            # > 1.0
            d = kernel_value_2D(k, 1.5, 0.5)
            @test d == 0.0
        end

        @testset "2D derivative" begin
            k = Quintic()
            # < 1/3
            d = kernel_deriv_2D(k, 0.3, 0.5)
            @test d ≈ -0.8389284295663888
            # < 2/3
            d = kernel_deriv_2D(k, 0.5, 0.5)
            @test d ≈ -95.94285222098028
            # < 1.0
            d = kernel_deriv_2D(k, 0.8, 0.5)
            @test d ≈ -0.010194573321313068
            # > 1.0
            d = kernel_deriv_2D(k, 1.5, 0.5)
            @test d == 0.0
        end

        @testset "2D bias correction" begin
            k = Quintic()
            d = bias_correction_2D(k, 1.0, 0.0, 0.0)
            @test d ≈  1.0  
        end

        @testset "3D value" begin
            k = Quintic()
            # < 1/3
            d = kernel_value_3D(k, 0.3, 0.5)
            @test d ≈ 0.2791208661307549
            # < 2/3    
            d = kernel_value_3D(k, 0.5, 0.5)
            @test d ≈ 0.06630419797168217
            # < 1.0
            d = kernel_value_3D(k, 0.8, 0.5)
            @test d ≈ 0.0006961437210839494
            # > 1.0
            d = kernel_value_3D(k, 1.5, 0.5)
            @test d == 0.0
        end

        @testset "3D derivative" begin
            k = Quintic()
            # < 1/3
            d = kernel_deriv_3D(k, 0.3, 0.5)
            @test d ≈ -0.7160853380941675
            # < 2/3
            d = kernel_deriv_3D(k, 0.5, 0.5)
            @test d ≈ -81.89407743147959
            # < 1.0
            d = kernel_deriv_3D(k, 0.8, 0.5)
            @test d ≈ -0.008701796513549367
            # > 1.0
            d = kernel_deriv_3D(k, 1.5, 0.5)
            @test d == 0.0
        end

        @testset "3D bias correction" begin
            k = Quintic()
            d = bias_correction_3D(k, 1.0, 0.0, 0.0)
            @test d ≈  1.0  
        end
    end

    @testset "Wendland C2" begin
        
        @testset "1D value" begin
            k = WendlandC2()
            # < 1.0
            d = kernel_value_1D(k, 0.5, 0.5)
            @test d ≈ 0.1953125
            # > 1.0
            d = kernel_value_1D(k, 1.5, 0.5)
            @test d == 0.0
        end

        @testset "1D derivative" begin
            k = WendlandC2()
            # < 1.0
            d = kernel_deriv_1D(k, 0.5, 0.5)
            @test d ≈ -0.9375
            # > 1.0
            d = kernel_deriv_1D(k, 1.5, 0.5)
            @test d == 0.0
        end

        @testset "1D bias correction" begin
            k = WendlandC2()
            d = bias_correction_1D(k, 1.0, 1.0, 0.5)

            @test d ≈ 0.99540625
        end

        @testset "2D value" begin
            k = WendlandC2()
            # < 1.0
            d = kernel_value_2D(k, 0.5, 0.5)
            @test d ≈ 0.10444543140405632
            # > 1.0
            d = kernel_value_2D(k, 1.5, 0.5)
            @test d == 0.0
        end

        @testset "2D derivative" begin
            k = WendlandC2()
            # < 1.0
            d = kernel_deriv_2D(k, 0.5, 0.5)
            @test d ≈ -0.3481514380135211
            # > 1.0
            d = kernel_deriv_2D(k, 1.5, 0.5)
            @test d == 0.0
        end

        @testset "2D bias correction" begin
            k = WendlandC2()
            d = bias_correction_2D(k, 1.0, 1.0, 0.5)

            @test d ≈ 0.991811478177922
        end

        @testset "3D value" begin
            k = WendlandC2()
            # < 1.0
            d = kernel_value_3D(k, 0.5, 0.5)
            @test d ≈ 0.07833407355304224
            # > 1.0
            d = kernel_value_3D(k, 1.5, 0.5)
            @test d == 0.0
        end

        @testset "3D derivative" begin
            k = WendlandC2()
            # < 1.0
            d = kernel_deriv_3D(k, 0.5, 0.5)
            @test d ≈ -0.26111357851014083
            # > 1.0
            d = kernel_deriv_3D(k, 1.5, 0.5)
            @test d == 0.0
        end

        @testset "3D bias correction" begin
            k = WendlandC2()
            d = bias_correction_3D(k, 1.0, 1.0, 0.5)

            @test d ≈ 0.987717217266883
        end
    end

    @testset "Wendland C4" begin
        
        @testset "1D value" begin
            k = WendlandC4()
            # < 1.0
            d = kernel_value_1D(k, 0.5, 0.5)
            @test d ≈ 0.12890625
            # > 1.0
            d = kernel_value_1D(k, 1.5, 0.5)
            @test d == 0.0
        end

        @testset "1D derivative" begin
            k = WendlandC4()
            # < 1.0
            d = kernel_deriv_1D(k, 0.5, 0.5)
            @test d ≈ -0.4921875
            # > 1.0
            d = kernel_deriv_1D(k, 1.5, 0.5)
            @test d == 0.0
        end

        @testset "1D bias correction" begin
            k = WendlandC4()
            d = bias_correction_1D(k, 1.0, 1.0, 0.5)

            @test d ≈ 0.9992541503078153
        end

        @testset "2D value" begin
            k = WendlandC4()
            # < 1.0
            d = kernel_value_2D(k, 0.5, 0.5)
            @test d ≈ 0.07740152505836316
            # > 1.0
            d = kernel_value_2D(k, 1.5, 0.5)
            @test d == 0.0
        end

        @testset "2D derivative" begin
            k = WendlandC4()
            # < 1.0
            d = kernel_deriv_2D(k, 0.5, 0.5)
            @test d ≈ -0.37301939787162974
            # > 1.0
            d = kernel_deriv_2D(k, 1.5, 0.5)
            @test d == 0.0
        end

        @testset "2D bias correction" begin
            k = WendlandC4()
            d = bias_correction_2D(k, 1.0, 1.0, 0.5)

            @test d ≈ 0.9985755320162227
        end

        @testset "3D value" begin
            k = WendlandC4()
            # < 1.0
            d = kernel_value_3D(k, 0.5, 0.5)
            @test d ≈ 0.06651693559703084
            # > 1.0
            d = kernel_value_3D(k, 1.5, 0.5)
            @test d == 0.0
        end

        @testset "3D derivative" begin
            k = WendlandC4()
            # < 1.0
            d = kernel_deriv_3D(k, 0.5, 0.5)
            @test d ≈ -0.3205635450459318
            # > 1.0
            d = kernel_deriv_3D(k, 1.5, 0.5)
            @test d == 0.0
        end

        @testset "3D bias correction" begin
            k = WendlandC4()
            d = bias_correction_3D(k, 1.0, 1.0, 0.5)

            @test d ≈ 0.9975516956528828
        end
    end

    @testset "Wendland C6" begin

        @testset "1D value" begin
            k = WendlandC6()
            # < 1.0
            d = kernel_value_1D(k, 0.5, 0.5)
            @test d ≈ 0.0797271728515625
            # > 1.0
            d = kernel_value_1D(k, 1.5, 0.5)
            @test d == 0.0
        end

        @testset "1D derivative" begin
            k = WendlandC6()
            # < 1.0
            d = kernel_deriv_1D(k, 0.5, 0.5)
            @test d ≈ -0.417938232421875
            # > 1.0
            d = kernel_deriv_1D(k, 1.5, 0.5)
            @test d == 0.0
        end

        @testset "1D bias correction" begin
            k = WendlandC6()
            d = bias_correction_1D(k, 1.0, 1.0, 0.5)

            @test d ≈ 0.9997781503817778
        end

        @testset "2D value" begin
            k = WendlandC6()
            # < 1.0
            d = kernel_value_2D(k, 0.5, 0.5)
            @test d ≈ 0.052822211162893276
            # > 1.0
            d = kernel_value_2D(k, 1.5, 0.5)
            @test d == 0.0
        end

        @testset "2D derivative" begin
            k = WendlandC6()
            # < 1.0
            d = kernel_deriv_2D(k, 0.5, 0.5)
            @test d ≈ -0.3238607700806899
            # > 1.0
            d = kernel_deriv_2D(k, 1.5, 0.5)
            @test d == 0.0
        end

        @testset "2D bias correction" begin
            k = WendlandC6()
            d = bias_correction_2D(k, 1.0, 1.0, 0.5)

            @test d ≈ 0.9995421822101074
        end

        @testset "3D value" begin
            k = WendlandC6()
            # < 1.0
            d = kernel_value_3D(k, 0.5, 0.5)
            @test d ≈ 0.05055250677698771
            # > 1.0
            d = kernel_value_3D(k, 1.5, 0.5)
            @test d == 0.0
        end

        @testset "3D derivative" begin
            k = WendlandC6()
            # < 1.0
            d = kernel_deriv_3D(k, 0.5, 0.5)
            @test d ≈ -0.30994487761628525
            # > 1.0
            d = kernel_deriv_3D(k, 1.5, 0.5)
            @test d == 0.0
        end

        @testset "3D bias correction" begin
            k = WendlandC6()
            d = bias_correction_3D(k, 1.0, 1.0, 0.5)

            @test d ≈ 0.9991237081365336
        end
    end

    @testset "Wendland C8" begin

        @testset "2D value" begin
            k = WendlandC8()
            # < 1.0
            d = kernel_value_2D(k, 0.5, 0.5)
            @test d ≈ 0.034310013366734275
            # > 1.0
            d = kernel_value_2D(k, 1.5, 0.5)
            @test d == 0.0
        end

        @testset "2D derivative" begin
            k = WendlandC8()
            # < 1.0
            d = kernel_deriv_2D(k, 0.5, 0.5)
            @test d ≈ -0.5179270722948169
            # > 1.0
            d = kernel_deriv_2D(k, 1.5, 0.5)
            @test d == 0.0
        end

        @testset "2D bias correction" begin
            k = WendlandC8()
            d = bias_correction_2D(k, 1.0, 1.0, 0.5)

            @test d ≈ 1.0
        end

        @testset "3D value" begin
            k = WendlandC8()
            # < 1.0
            d = kernel_value_3D(k, 0.5, 0.5)
            @test d ≈ 0.034310013366734275
            # > 1.0
            d = kernel_value_3D(k, 1.5, 0.5)
            @test d == 0.0
        end

        @testset "3D derivative" begin
            k = WendlandC8()
            # < 1.0
            d = kernel_deriv_3D(k, 0.5, 0.5)
            @test d ≈ -0.5179270722948169
            # > 1.0
            d = kernel_deriv_3D(k, 1.5, 0.5)
            @test d == 0.0
        end

        @testset "3D bias correction" begin
            k = WendlandC8()
            d = bias_correction_3D(k, 1.0, 1.0, 0.5)

            @test d ≈ 1.0
        end
    end

    @testset "Multiple Dispatch" begin

        @testset "kernel value" begin
            k = WendlandC6()
            @test 𝒲₁(k, 0.5, 1.0) ≈ kernel_value_1D(k, 0.5, 1.0)
            @test 𝒲₂(k, 0.5, 1.0) ≈ kernel_value_2D(k, 0.5, 1.0)
            @test 𝒲₃(k, 0.5, 1.0) ≈ kernel_value_3D(k, 0.5, 1.0)
        end

        @testset "kernel derivative" begin
            k = WendlandC6()
            @test ∇𝒲₁(k, 0.5, 1.0) ≈ kernel_deriv_1D(k, 0.5, 1.0)
            @test ∇𝒲₂(k, 0.5, 1.0) ≈ kernel_deriv_2D(k, 0.5, 1.0)
            @test ∇𝒲₃(k, 0.5, 1.0) ≈ kernel_deriv_3D(k, 0.5, 1.0)
        end

        @testset "bias correction" begin
            k = WendlandC6()
            @test δρ₁(k, 1.0, 1.0, 0.5) ≈ bias_correction_1D(k, 1.0, 1.0, 0.5)
            @test δρ₂(k, 1.0, 1.0, 0.5) ≈ bias_correction_2D(k, 1.0, 1.0, 0.5)
            @test δρ₃(k, 1.0, 1.0, 0.5) ≈ bias_correction_3D(k, 1.0, 1.0, 0.5)
        end
        
    end

end