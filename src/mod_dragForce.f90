!
!   Module to estimate drag force on a particle in Stokes flow with respect to neighbouring particles
!
MODULE mDragModel
     
    use mNeuralNets
    use mCommon
    use mRotation
    IMPLICIT NONE
!
!    Parameters
!
! ----------------------------------------------------------------------
!

!
!    Types
!

!
! ----------------------------------------------------------------------
!
!
!   Variables
!

    CONTAINS

    subroutine testDragForce(streamwiseFirst, lateralFirst, streamwiseSecond, lateralSecond)
        real(rk) targetParticle(3)
        real(rk) p1(3),p2(3),p3(3),p4(3),p5(3)
        real(rk) targetParticleVelocity(3)
        real(rk) flowVelocityAtTarget(3)
        real(rk) viscosity, diameter
        real(rk) dragForce(3)
        TYPE (NeuralNet) :: streamwiseFirst, lateralFirst
        TYPE (NeuralNet) :: streamwiseSecond, lateralSecond 

    
        targetParticle =         (/ 0.023_rk,  0.0051_rk,  0.05_rk /)
        p1 =                     (/ 0.022_rk,  0.00066_rk, 0.045_rk /)
        p2 =                     (/ 0.027_rk,  0.0055_rk,  0.041_rk /)
        p3 =                     (/ 0.033_rk,  0.0055_rk,  0.047_rk /)
        p4 =                     (/ 0.020_rk,  0.0037_rk,  0.034_rk /)
        p5 =                     (/ 0.010_rk,  0.0043_rk,  0.039_rk /)    
        targetParticleVelocity = (/ 0.0005_rk, 0.0001_rk,  -0.0002_rk /)
        flowVelocityAtTarget =   (/ 0.0001_rk,  -0.0004_rk,  0.0003_rk /)
        viscosity = 0.01_rk
        diameter = 1e-4_rk
    
        call getForce(targetParticle,p1,p2,p3,p4,p5,targetParticleVelocity,flowVelocityAtTarget,diameter,viscosity, dragForce,& 
                      streamwiseFirst, lateralFirst, streamwiseSecond, lateralSecond)

   
    end subroutine    


    subroutine getForce(targetParticle,p1,p2,p3,p4,p5,targetParticleVelocity,flowVelocityAtTarget,diameter,viscosity,dragForce,&
                        streamwiseFirst, lateralFirst, streamwiseSecond, lateralSecond)
 
        TYPE (NeuralNet) :: streamwiseFirst, lateralFirst
        TYPE (NeuralNet) :: streamwiseSecond, lateralSecond 
        real(rk) targetParticle(3)
        real(rk) p1(3),p2(3),p3(3),p4(3),p5(3)
        real(rk) targetParticleVelocity(3)
        real(rk) flowVelocityAtTarget(3)
        real(rk) dragForce(3),beta(3)
        real(rk) diameter, viscosity
        real(rk) relativeVelocity(3)

        real(rk) p1sph(3),p2sph(3),p3sph(3),p4sph(3),p5sph(3)
        real(rk) rotationMatrix(3,3)
        integer i,j
        !
        !  Get rotation matrix, so velocity is in y directon
        !
        call constructRotationMatrix(targetParticleVelocity,flowVelocityAtTarget,rotationMatrix)
        if (debug) then
            print *,"rotationMatrix"
            print '(3F10.4)', ((rotationMatrix(i, j), j=1, 3), i=1, 3)
        end if

        !
        !  Transform particle coordinates to spherical cs
        !
        call getSPH(rotationMatrix,targetParticle,p1,p2,p3,p4,p5,p1sph,p2sph,p3sph,p4sph,p5sph,diameter)
        !
        !  Use neural network to get Stokes force multiplyer beta
        !
        call getBeta(p1sph,p2sph,p3sph,p4sph,p5sph,beta, &
                     streamwiseFirst, lateralFirst, streamwiseSecond, lateralSecond)
        if (debug) print *,"beta",beta
        !
        !  Transform relative velocity to local cs
        !

        relativeVelocity = flowVelocityAtTarget - targetParticleVelocity
        if (debug) print *,"rVel GCS",relativeVelocity        
        relativeVelocity = matmul(rotationMatrix,relativeVelocity)
        if (debug) print *,"rVel LCS",relativeVelocity        
        !
        !  Calculate drag force in local cs
        !
        dragForce = beta * 3.0_rk * pi * viscosity * diameter * relativeVelocity(2) 
        if (debug) print *,"dragForce LCS",dragForce
        !
        ! Rotate back to global cs
        !
        dragForce = matmul(transpose(rotationMatrix),dragForce)
        if (debug) print *,"dragForce GCS",dragForce


    end subroutine


    subroutine getSPH(rotationMatrix,targetParticle,p1,p2,p3,p4,p5,p1sph,p2sph,p3sph,p4sph,p5sph,diameter)

        real(rk) p1(3),p2(3),p3(3),p4(3),p5(3)
        real(rk) p1lcs(3),p2lcs(3),p3lcs(3),p4lcs(3),p5lcs(3)
        real(rk) p1sph(3),p2sph(3),p3sph(3),p4sph(3),p5sph(3)
        real(rk) rotationMatrix(3,3),diameter
        real(rk) targetParticle(3)
        real(rk) d1,d2,d3,d4,d5

        d1 = distance(p1,targetParticle)
        d2 = distance(p2,targetParticle)
        d3 = distance(p3,targetParticle)
        d4 = distance(p4,targetParticle)
        d5 = distance(p5,targetParticle)

        if (debug) then
            print *, "Particles x,y,z & distance"
            print *, p1(1), p1(2), p1(3), d1
            print *, p2(1), p2(2), p2(3), d2
            print *, p3(1), p3(2), p3(3), d3
            print *, p4(1), p4(2), p4(3), d4
            print *, p5(1), p5(2), p5(3), d5
        end if

        !
        ! Relative coordinate system with target particle at (0,0,0)
        !
        p1lcs = p1 - targetParticle
        p2lcs = p2 - targetParticle
        p3lcs = p3 - targetParticle
        p4lcs = p4 - targetParticle
        p5lcs = p5 - targetParticle
        !
        ! Rotate
        !
        p1lcs = matmul(rotationMatrix,p1lcs)
        p2lcs = matmul(rotationMatrix,p2lcs)
        p3lcs = matmul(rotationMatrix,p3lcs)        
        p4lcs = matmul(rotationMatrix,p4lcs)
        p5lcs = matmul(rotationMatrix,p5lcs)
        if (debug) then
            print *, "Rotated particles lcs"
            print *, p1lcs(1), p1lcs(2), p1lcs(3)
            print *, p2lcs(1), p2lcs(2), p2lcs(3)
            print *, p3lcs(1), p3lcs(2), p3lcs(3)
            print *, p4lcs(1), p4lcs(2), p4lcs(3)
            print *, p5lcs(1), p5lcs(2), p5lcs(3)
        end if        
        !
        ! Spherical coordiante system
        !
        p1sph(1) = d1 / diameter ! r
        p2sph(1) = d2 / diameter
        p3sph(1) = d3 / diameter
        p4sph(1) = d4 / diameter
        p5sph(1) = d5 / diameter

        p1sph(2) = acos( p1lcs(3) / d1 ) ! theta
        p2sph(2) = acos( p2lcs(3) / d2 ) ! theta
        p3sph(2) = acos( p3lcs(3) / d3 ) ! theta
        p4sph(2) = acos( p4lcs(3) / d4 ) ! theta
        p5sph(2) = acos( p5lcs(3) / d5 ) ! theta

        p1sph(3) = atan2( p1lcs(2) , p1lcs(1) ) ! phi
        p2sph(3) = atan2( p2lcs(2) , p2lcs(1) ) ! phi
        p3sph(3) = atan2( p3lcs(2) , p3lcs(1) ) ! phi
        p4sph(3) = atan2( p4lcs(2) , p4lcs(1) ) ! phi
        p5sph(3) = atan2( p5lcs(2) , p5lcs(1) ) ! phi
        if (debug) then
            print *, "Spherical cs particles: r,theta,phi"
            print *, p1sph(1), p1sph(2), p1sph(3)
            print *, p2sph(1), p2sph(2), p2sph(3)
            print *, p3sph(1), p3sph(2), p3sph(3)
            print *, p4sph(1), p4sph(2), p4sph(3)
            print *, p5sph(1), p5sph(2), p5sph(3)
            print *,""
        end if  

    end subroutine


    subroutine getBeta(p1sph,p2sph,p3sph,p4sph,p5sph,beta, &
                      streamwiseFirst, lateralFirst, streamwiseSecond, lateralSecond)

        TYPE (NeuralNet) :: streamwiseFirst, lateralFirst
        TYPE (NeuralNet) :: streamwiseSecond, lateralSecond

        real(rk) p1sph(3),p2sph(3),p3sph(3),p4sph(3),p5sph(3)
        real(rk) beta(3)
        real(rk), allocatable :: inputLayer(:),outputLayer(:)

        beta = (/ 0.0_rk, 0.0_rk, 0.0_rk /) ! result = % of stokes drag
        !
        !  Use the streamwiseFirst neural net for p1 & p2
        !
        allocate(outputLayer(streamwiseFirst%outLength))
        allocate(inputLayer(streamwiseFirst%inpLength))

        inputLayer(1) = p1sph(1)
        inputLayer(2) = p1sph(2)
        inputLayer(3) = p1sph(3)
        inputLayer(4) = p2sph(1)
        inputLayer(5) = p2sph(2)
        inputLayer(6) = p2sph(3)
        
        call nnEval(streamwiseFirst,inputLayer,outputLayer)
        beta(2) = beta(2) + outputLayer(1)  / 24.0_rk 
        deallocate(inputLayer,outputLayer)
        !
        !  Use the streamwiseSecond neural net for p3
        !
        allocate(outputLayer(streamwiseSecond%outLength))
        allocate(inputLayer(streamwiseSecond%inpLength))

        inputLayer(1) = p3sph(1)
        inputLayer(2) = p3sph(2)
        inputLayer(3) = p3sph(3)
        
        call nnEval(streamwiseSecond,inputLayer,outputLayer)
        beta(2) = beta(2) + outputLayer(1)  / 24.0_rk
        deallocate(inputLayer,outputLayer)        
        !
        !  Use the streamwiseSecond neural net for p4
        !
        allocate(outputLayer(streamwiseSecond%outLength))
        allocate(inputLayer(streamwiseSecond%inpLength))

        inputLayer(1) = p4sph(1)
        inputLayer(2) = p4sph(2)
        inputLayer(3) = p4sph(3)
        
        call nnEval(streamwiseSecond,inputLayer,outputLayer)
        beta(2) = beta(2) + outputLayer(1)  / 24.0_rk
        deallocate(inputLayer,outputLayer)    
        !
        !  Use the streamwiseSecond neural net for p3
        !
        allocate(outputLayer(streamwiseSecond%outLength))
        allocate(inputLayer(streamwiseSecond%inpLength))

        inputLayer(1) = p5sph(1)
        inputLayer(2) = p5sph(2)
        inputLayer(3) = p5sph(3)
        
        call nnEval(streamwiseSecond,inputLayer,outputLayer)
        beta(2) = beta(2) + outputLayer(1)  / 24.0_rk
        deallocate(inputLayer,outputLayer)    
        !
        ! Lateral
        !

        !
        !  Use the streamwiseFirst neural net for p1 & p2
        !
        allocate(outputLayer(lateralFirst%outLength))
        allocate(inputLayer(lateralFirst%inpLength))

        inputLayer(1) = p1sph(1)
        inputLayer(2) = p1sph(2)
        inputLayer(3) = p1sph(3)
        inputLayer(4) = p2sph(1)
        inputLayer(5) = p2sph(2)
        inputLayer(6) = p2sph(3)
        
        call nnEval(lateralFirst,inputLayer,outputLayer)

        beta(1) = beta(1) + outputLayer(1)  / 24.0_rk 
        beta(3) = beta(3) + outputLayer(2)  / 24.0_rk 
        deallocate(inputLayer,outputLayer)
        !
        !  Use the lateralSecond neural net for p3
        !
        allocate(outputLayer(lateralSecond%outLength))
        allocate(inputLayer(lateralSecond%inpLength))

        inputLayer(1) = p3sph(1)
        inputLayer(2) = p3sph(2)
        inputLayer(3) = p3sph(3)
        
        call nnEval(lateralSecond,inputLayer,outputLayer)

        beta(1) = beta(1) + outputLayer(1)  / 24.0_rk
        beta(3) = beta(3) + outputLayer(2)  / 24.0_rk
        deallocate(inputLayer,outputLayer)        
        !
        !  Use the lateralSecond neural net for p4
        !
        allocate(outputLayer(lateralSecond%outLength))
        allocate(inputLayer(lateralSecond%inpLength))

        inputLayer(1) = p4sph(1)
        inputLayer(2) = p4sph(2)
        inputLayer(3) = p4sph(3)
        
        call nnEval(lateralSecond,inputLayer,outputLayer)
        beta(1) = beta(1) + outputLayer(1)  / 24.0_rk
        beta(3) = beta(3) + outputLayer(2)  / 24.0_rk
        deallocate(inputLayer,outputLayer)    
        !
        !  Use the lateralSecond neural net for p5
        !
        allocate(outputLayer(lateralSecond%outLength))
        allocate(inputLayer(lateralSecond%inpLength))

        inputLayer(1) = p5sph(1)
        inputLayer(2) = p5sph(2)
        inputLayer(3) = p5sph(3)
        
        call nnEval(lateralSecond,inputLayer,outputLayer)
        beta(1) = beta(1) + outputLayer(1)  / 24.0_rk
        beta(3) = beta(3) + outputLayer(2)  / 24.0_rk
        deallocate(inputLayer,outputLayer)   

    end

    function distance(n1, n2) result(d)
        real(rk), intent(in) :: n1(3), n2(3)
        real(rk) :: d

        d = sqrt( ( n1(1) - n2(1) )**2 + (n1(2) - n2(2))**2 + (n1(3) - n2(3))**2 )

    end function


    subroutine constructRotationMatrix(targetParticleVelocity,flowVelocityAtTarget,rotationMatrix)
    
        real(rk) targetParticleVelocity(3)
        real(rk) flowVelocityAtTarget(3)
        real(rk) rotationMatrix(3,3)
        real(rk) relativeVelocity(3)
        real(rk) relativeVelocityMag
        real(rk) e2(3),e2prime(3)
        
        relativeVelocity = flowVelocityAtTarget - targetParticleVelocity
        relativeVelocityMag = sqrt ( relativeVelocity(1)**2 + relativeVelocity(2)**2 + relativeVelocity(3)**2 )
        e2(1) = 0.0_rk
        e2(2) = 1.0_rk
        e2(3) = 0.0_rk
        if ( relativeVelocityMag .GT. vSmall ) then
            e2prime = relativeVelocity / relativeVelocityMag
            rotationMatrix = rotationTensor(e2Prime, e2)
        else 
            rotationMatrix(1,1) = 1.0_rk
            rotationMatrix(1,2) = 0.0_rk
            rotationMatrix(1,3) = 0.0_rk
            rotationMatrix(2,1) = 0.0_rk
            rotationMatrix(2,2) = 1.0_rk
            rotationMatrix(2,3) = 0.0_rk
            rotationMatrix(3,1) = 0.0_rk
            rotationMatrix(3,2) = 0.0_rk
            rotationMatrix(3,3) = 1.0_rk
        end if
    
    end subroutine


END MODULE mDragModel