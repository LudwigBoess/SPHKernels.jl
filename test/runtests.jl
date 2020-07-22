using SPHKernels, Test

@testset "SPH Kernels" begin

    @testset "Cubic Spline" begin

        @testset "2D value" begin
            k = Cubic()
            # < 0.5
            d = kernel_value_2D(k, 0.4, 0.5)
            @test d ≈ 0.26992678348385446       
            # < 1.0
            d = kernel_value_2D(k, 0.5, 0.5)
            @test d ≈ 0.15915494309189535
            # > 1.0
            d = kernel_value_2D(k, 1.5, 0.5)
            @test d == 0.0
        end

        @testset "2D derivative" begin
            k = Cubic()
            # < 0.5
            d = kernel_deriv_2D(k, 0.4, 0.5)
            @test d ≈ -0.6111549814728781      
            # < 1.0
            d = kernel_deriv_2D(k, 0.5, 0.5)
            @test d ≈ -0.477464829275686
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
            @test d ≈ 0.5016563806256541
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
            @test d ≈-95.94285222098028
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

    @testset "Wendland C4" begin
        
        @testset "2D value" begin
            k = WendlandC4()
            # < 1.0
            d = kernel_value_2D(k, 0.5, 0.5)
            @test d ≈ 0.15480305011672632
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
            @test d ≈ 0.13303387119406168
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

end