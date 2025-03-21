subroutine randomlyInsertParticles(particles, simulationDomainSideLength, n)
    use mCommon
    implicit none
    integer n
    real(rk), intent(inout) :: particles(n,3)
    real(rk), intent(in) :: simulationDomainSideLength
    real(rk) rnd
    integer :: i
    call random_seed()
   
    do i = 1, size(particles, 1)
        call random_number(rnd); particles(i, 1) = rnd * simulationDomainSideLength
        call random_number(rnd); particles(i, 2) = rnd * simulationDomainSideLength
        call random_number(rnd); particles(i, 3) = rnd * simulationDomainSideLength
        !particles(i, :) = particles(i, :) * simulationDomainSideLength
    end do
end subroutine randomlyInsertParticles




subroutine find5closest(particles,p1,p2,p3,p4,p5,targetParticle,id,n,particleVolume,d5)
    use mCommon
    implicit none
    integer i,id,n
    real(rk) p1(3),p2(3),p3(3),p4(3),p5(3)
    real(rk) particles(n,3)
    real(rk) targetParticle(3),p(3),d5
    real(rk), allocatable :: d2(:)
    real(rk) minDist,dist2,sphereRadius,particleVolume
    integer found

    targetParticle = particles(id,:)

    allocate(d2(n))
    do i=1,n
        p = particles(i,:)
        d2(i) = dist2(targetParticle,p)
    end do

    minDist = 1.0E15_rk
    do i=1,n
        if (i.ne.id) then
            if ( d2(i) .LT. minDist ) then
                minDist = d2(i)
                found = i
            end if
        end if
    end do
    p1 = particles(found,:)
    minDist = 1.0E15_rk
    d2(found) = minDist

    do i=1,n
        if (i.ne.id) then
            if ( d2(i) .LT. minDist ) then
                minDist = d2(i)
                found = i
            end if
        end if
    end do
    p2 = particles(found,:)
    minDist = 1.0E15_rk
    d2(found) = minDist

    do i=1,n
        if (i.ne.id) then
            if ( d2(i) .LT. minDist ) then
                minDist = d2(i)
                found = i
            end if
        end if
    end do
    p3 = particles(found,:)
    minDist = 1.0E15_rk
    d2(found) = minDist    

    do i=1,n
        if (i.ne.id) then
            if ( d2(i) .LT. minDist ) then
                minDist = d2(i)
                found = i
            end if
        end if
    end do
    p4 = particles(found,:)
    minDist = 1.0E15_rk
    d2(found) = minDist    

    do i=1,n
        if (i.ne.id) then
            if ( d2(i) .LT. minDist ) then
                minDist = d2(i)
                found = i
            end if
        end if
    end do
    p5 = particles(found,:)
    minDist = 1.0E15_rk
    d2(found) = minDist        

    ! 
    ! Rembemer the distance to 5th particle
    !
    d5 = sqrt(dist2(targetParticle,p5))
 
    deallocate(d2)
 
end subroutine


real(rk) function dist2(n1, n2) result(d)
    use mCommon
    real(rk), intent(in) :: n1(3), n2(3)

    d =  ( n1(1) - n2(1) )**2 + (n1(2) - n2(2))**2 + (n1(3) - n2(3))**2 

end function