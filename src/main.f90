program settling

    use mCommon
    use mDragModel

    real(rk), allocatable :: particles(:,:),particlesVelocity(:,:), d5(:)
    real(rk) :: volumeFraction, simulationDomainSideLength
    integer :: n,m,i,j,v,vv
    real(rk) t,dt,eps

    real(rk) flowVelocityAtTarget(3),particleCharacteristicTime
    real(rk) :: variance, stddev,averageVelocity,oldAverageVelocity,averageLocalVF

    real(rk) dragForce(3),gravity(3),bouyancy(3),totalForce(3),sumTotalForce(3),sumDragForce(3)
    real(rk) particleTerminalVelocity
    real(rk) diameter, viscosity, fluidDensity, particleDensity, gravConstant,particleMass,particleVolume  
    real(rk) p1(3),p2(3),p3(3),p4(3),p5(3)
    real(rk) targetParticle(3),targetParticleVelocity(3)
    character(255) vf,fname
    integer exPart13,exPart14,exPart15
    logical done

    TYPE (NeuralNet) :: streamwiseFirst, lateralFirst
    TYPE (NeuralNet) :: streamwiseSecond, lateralSecond
    !
    ! Get process ID
    !    
    !character(255) pid
    !write(pid, '(I0)') getpid()
    !call system("echo newTimeStep >> memory_usage.log")
    !call system("ps -o rss= -p " // trim(pid) // " >> memory_usage.log")    

    !
    ! Read in neural networks
    !
    streamwiseFirst = nnCreate("streamwiseFirst.nn")
    streamwiseSecond = nnCreate("streamwiseSecond.nn")

    lateralFirst = nnCreate("lateralFirst.nn")
    lateralSecond = nnCreate("lateralSecond.nn")
    !


!    call testDragForce(streamwiseFirst, lateralFirst, streamwiseSecond, lateralSecond)
!    stop

    flowVelocityAtTarget = (/ 0.0_rk, 0.0_rk, 0.0_rk /)
    diameter = 25E-6 !m
    viscosity = 1.81E-5_rk ! Pa.s (viscosity of air at 20C)
    fluidDensity = 1.225_rk ! kg/m^3 (density of air at 20C)
    particleDensity = 2500_rk ! kg/m^3
    gravConstant = 9.81 ! m/s2
    particleVolume=(4.0_rk / 3.0_rk) * pi * (diameter / 2.0_rk)**3
    particleMass = particleDensity *  particleVolume

    particleTerminalVelocity = 2.0_rk/9.0_rk * (diameter / 2.0_rk)**2 *gravConstant /(viscosity)*(particleDensity - fluidDensity)
    particleCharacteristicTime = particleMass /(3.0_rk * pi * viscosity * diameter)
    dt = particleCharacteristicTime/100.0_rk ! time step

    exPart13 = 33
    WRITE(fname,'(A,I0,A)') "singleParticle",exPart13,".dat"
    open(13,FILE=trim(fname),STATUS="UNKNOWN")
    WRITE(13,*) "j t dFx dfY dfZ tFx tFy tFz x y z vx vy vz vx/s vy/s vz/s lVf"

    exPart14 = 66
    WRITE(fname,'(A,I0,A)') "singleParticle",exPart14,".dat"
    open(14,FILE=trim(fname),STATUS="UNKNOWN")
    WRITE(14,*) "j t dFx dfY dfZ tFx tFy tFz x y z vx vy vz vx/s vy/s vz/s lVf"

    exPart15 = 99
    WRITE(fname,'(A,I0,A)') "singleParticle",exPart15,".dat"
    open(15,FILE=trim(fname),STATUS="UNKNOWN")
    WRITE(15,*) "j t dFx dfY dfZ tFx tFy tFz x y z vx vy vz vx/s vy/s vz/s lVf"    


    open(12,FILE="vFaV.dat",STATUS="UNKNOWN")

    write(12,*)  "flowVelocityAtTarget:", flowVelocityAtTarget
    write(12,*)  "diameter:", diameter
    write(12,*)  "viscosity:", viscosity
    write(12,*)  "fluidDensity:", fluidDensity
    write(12,*)  "particleDensity:", particleDensity
    write(12,*)  "gravConstant:", gravConstant
    write(12,*)  "particleVolume:", particleVolume
    write(12,*)  "particleMass:", particleMass
    write(12,*)  "particleTerminalVelocity:", particleTerminalVelocity
    write(12,*)  "particleCharacteristicTime:", particleCharacteristicTime
    write(12,*)  "dt:", dt
    write(12,*) &
    "j t lVf av/s stdev/s avVel stdev avTFx avTFy avTFz avTFx/g avTFy/g avTFz/g avDFx avDFy avDFz avDFx/g avDFy/g avDFz/g"

    eps = 1.0E-10

    vv = 5
    call get_command_argument(1, vf)
    if (len_trim(vf) .ne. 0) then
        read(vf, *) vv
    end if

    do v = vv,vv


        volumeFraction = 1.0_rk / 10.0**v
        write(12,*)  "volumeFraction:", volumeFraction
    
        n = 100000 ! number of particles
        m = 10000 ! number of time steps
        t = 0.0_rk ! time   
        oldAverageVelocity = 0.0_rk
        

        allocate(particles(n,3),particlesVelocity(n,3),d5(n))
        particlesVelocity = 0.0_rk

        simulationDomainSideLength = (n * particleVolume / volumeFraction)**(1.0_rk / 3.0_rk)
        call randomlyInsertParticles(particles, simulationDomainSideLength, n)
        
        done = .FALSE.
        j=0
        do while ( (done.EQV..FALSE.) .AND. (j.LT.m) ) 
            j=j+1
            !
            ! Loop over particles
            !
            sumTotalForce = 0.0_rk
            sumDragForce = 0.0_rk
            do i=1,n
                !
                ! Find closest 5 particles
                !
                call find5closest(particles,p1,p2,p3,p4,p5,targetParticle,i,n,particleVolume,d5(i))
                !
                ! Calculate forces
                !
                targetParticleVelocity = particlesVelocity(i,:)
                call getForce(targetParticle,p1,p2,p3,p4,p5,targetParticleVelocity,flowVelocityAtTarget, &
                diameter,viscosity,dragForce,streamwiseFirst, lateralFirst, streamwiseSecond, lateralSecond)

                ! Stokes version
                !dragForce(1) = 0.0_rk
                !dragForce(2) = 3.0_rk * pi * viscosity * diameter * ( flowVelocityAtTarget(2) - particlesVelocity(i,2) )
                !dragForce(3) = 0.0_rk

                gravity = (/ 0.0_rk, - gravConstant * particleMass, 0.0_rk /)
                bouyancy = (/ 0.0_rk, fluidDensity * gravConstant * particleVolume, 0.0_rk /)

                totalForce = dragForce + gravity + bouyancy
                sumTotalForce = sumTotalForce + totalForce
                sumDragForce = sumDragForce + dragForce
                !
                ! Move particles
                !
                particlesVelocity(i,:) = particlesVelocity(i,:) + dt * totalForce / particleMass
                particles(i,:) = particles(i,:) + particlesVelocity(i,:) * dt   
                !
                ! Apply boundary conditions
                !
                if (particles(i,1).GT.simulationDomainSideLength) particles(i,1)=particles(i,1)-simulationDomainSideLength
                if (particles(i,2).GT.simulationDomainSideLength) particles(i,2)=particles(i,2)-simulationDomainSideLength
                if (particles(i,3).GT.simulationDomainSideLength) particles(i,3)=particles(i,3)-simulationDomainSideLength

                if (particles(i,1).LT.0.0) particles(i,1)=particles(i,1)+simulationDomainSideLength
                if (particles(i,2).LT.0.0) particles(i,2)=particles(i,2)+simulationDomainSideLength
                if (particles(i,3).LT.0.0) particles(i,3)=particles(i,3)+simulationDomainSideLength

                !
                ! Export a single particle
                !
                if (i.EQ.exPart13) then                   
                    WRITE(13,*) j,t,dragForce, totalForce, particles(i,:),particlesVelocity(i,:), &
                                particlesVelocity(i,:)/particleTerminalVelocity, &
                                15.0_rk * particleVolume / ( 4.0_rk * pi * ( d5(i) ) **3 )
                end if
                if (i.EQ.exPart14) then                   
                    WRITE(14,*) j,t,dragForce, totalForce, particles(i,:),particlesVelocity(i,:), &
                                particlesVelocity(i,:)/particleTerminalVelocity, &
                                15.0_rk * particleVolume / ( 4.0_rk * pi * ( d5(i) ) **3 )
                end if
                if (i.EQ.exPart15) then                   
                    WRITE(15,*) j,t,dragForce, totalForce, particles(i,:),particlesVelocity(i,:), &
                                particlesVelocity(i,:)/particleTerminalVelocity, &
                                15.0_rk * particleVolume / ( 4.0_rk * pi * ( d5(i) ) **3 )
                end if




            end do
            !
            ! Advance time
            !
            t = t + dt


            !
            ! Get average volume fraction
            !
            averageLocalVF = 0.0_rk
            do i=1,n
                averageLocalVF = averageLocalVF + d5(i)            
            end do
            averageLocalVF = averageLocalVF / dble(n)
            averageLocalVF = 15.0_rk * particleVolume / ( 4.0_rk * pi * ( averageLocalVF ) **3 )

            !
            ! Get average velocity
            !
            oldAverageVelocity = averageVelocity
            averageVelocity = 0.0_rk
            do i=1,n
                averageVelocity = averageVelocity + particlesVelocity(i,2)            
            end do
            averageVelocity = averageVelocity / dble(n)

            ! Calculate standard deviation
            variance = 0.0_rk
            do i=1,n
                variance = variance + (particlesVelocity(i,2) - averageVelocity)**2
            end do
            variance = variance / dble(n)
            stddev = sqrt(variance)

            if (abs(oldAverageVelocity-averageVelocity).LT.eps) done=.TRUE.
 
            write(12,*) j, t, &
                        averageLocalVF, &
                        abs(averageVelocity)/particleTerminalVelocity, &
                        stddev/particleTerminalVelocity, &
                        abs(averageVelocity),stddev, &
                        sumTotalForce/(dble(n)), sumTotalForce/(dble(n) * (gravity(2) + bouyancy(2)) ), &
                        sumDragForce/(dble(n)), sumDragForce/(dble(n) * (gravity(2) + bouyancy(2)) )
        end do

        deallocate(particles,particlesVelocity)    
    end do
    close(12)
    close(13)
    close(14)
    close(15)

end

