module Solvers

    contains

        subroutine SolveLinear(A,b,x,Solvable,argVerbose)
        
            complex,dimension(1:,1:) :: A
        
            complex,dimension(1:) :: b

            complex,dimension(1:size(b)) :: x

            complex,dimension(1:size(A,1),1:(size(A,1)+1)) :: Ab

            logical,optional :: Solvable,argVerbose

            logical :: swapFound,Verbose

            if(present(Solvable))then
                Solvable=.true.
            endif

            if(present(argVerbose))then
                Verbose=argVerbose
            else
                Verbose=.false.
            endif

            !Populate the Ab matrix with the matrix A and vector b
            do i=1,size(b)
                Ab(i,size(A,1)+1)=b(i)
            enddo
            do i=1,size(A,1)
                do j=1,size(A,2)
                    Ab(i,j)=A(i,j)
                enddo
            enddo

            !If any elements are too small, enlarge them all until proper numbers are reached
            ! do i=1,size(Ab,1)
                ! do j=1,size(Ab,2)
                    ! if(abs(Ab(i,j)).lt.1)then
                        ! do while(abs(Ab(i,j)).lt.1)
                            ! Ab=Ab*2
                        ! enddo
                    ! endif
                ! enddo
            ! enddo

            if(Verbose)then
                do i=1,size(Ab,1)
                    print *, Ab(i,:)
                enddo
            endif

            !Main oerations
            do i=1,size(Ab,1)
                !Handle zeros on diagonal
                if(Ab(i,i).eq.(0.,0.))then
                    if(Verbose) print *, 'Element',i,i,'is 0'
                    !Assume no swap found
                    swapFound=.false.
                    if(i.lt.size(Ab,1))then
                        do j=i,size(Ab,1)
                            if(Ab(j,i).ne.(0.,0.))then
                                swapFound=.true.
                                !Swapping
                                swappingBlock : block
                                    complex,dimension(1:size(Ab,2)) :: tempAb
                                    tempAb=Ab(i,:)
                                    Ab(i,:)=Ab(j,:)
                                    Ab(j,:)=tempAb
                                end block swappingBlock
                            endif
                        enddo
                    endif
                !If no swap found, skip the element
                    if(.not.swapFound)then
                        if(Verbose)then
                            print *, 'Module Solvers subroutine SolveLinear: zero element'
                            print *, Solvable
                        endif
                        cycle
                    endif
                endif
                !Gaussian elimination
                !Careful to do from furthest element so as not to make the diagonal 1 before elimination
                do j=size(Ab,2),1,-1
                    !Divide by the diagonal elements all elements in the same row
                    Ab(i,j) = Ab(i,j)/Ab(i,i)
                    if(Verbose)then
                        do k=1,size(Ab,1)
                            print *, Ab(k,:)
                        enddo
                    endif
                enddo
                !ELiminate column
                do j=1,size(Ab,1)
                    if(j.ne.i)then
                        !Do from furthest element first
                        do k=size(Ab,2),1,-1
                            Ab(j,k) = Ab(j,k)-Ab(i,k)*Ab(j,i)
                            if(Verbose)then
                                do l=1,size(Ab,1)
                                    print *, Ab(l,:)
                                enddo
                            endif
                        enddo
                    endif
                enddo
            enddo

            !After elimination complete check the zero element
            !Make sure to check from bottom up
            do i=size(A,1),1,-1
                if(Ab(i,i).eq.(0.,0.))then
                    if(Ab(i,size(Ab,2)).ne.(0.,0.))then
                        if(present(Solvable))then
                            Solvable=.false.
                            if(Verbose)then
                                print *, "Module Solvers subroutine SolveLinear: Not Solvable!"
                            endif
                        endif
                        return
                    else
                        cycle
                    endif
                endif
            enddo

            do i=1,size(Ab,1)
                print *, Ab(i,:)
            enddo

            !The solution is on the last column
            x=Ab(:,size(Ab,2))
        
        end subroutine SolveLinear

end module Solvers
