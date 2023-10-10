module Functions

    contains

        function KroneckerDelta(i,j)result(delta)

            integer :: i,j,delta

            if(i.eq.j)then
                delta=1
                return
            else
                delta=0
                return
            endif

        end function KroneckerDelta

        function LeviCivita(Indices)result(e)

            !Indices must be consecutive natural number starting from 1
            integer,dimension(1:) :: Indices

            integer :: e

            integer :: swapCount

            integer :: sumIndices

            integer :: twiceSum1toN

            swapCount=0
            sumIndices=0
            twiceSum1toN=0

            !This is a check for duplicate indices
            do i=1,size(Indices)
                !Convert indices into consecutive natural number from 1 via modulo function
                Indices(i)=modulo(Indices(i),size(Indices)+1)
                !Get the sum of indices
                sumIndices=sumIndices+Indices(i)
            enddo
            ! print *, sumIndices

            !Compare to the sum of consecutive natural number from 1
            twiceSum1toN=size(Indices)*(size(Indices)+1)
            ! print *, twiceSum1toN

            !If not equal return 0
            if((2*sumIndices).ne.twiceSum1toN)then
                e=0
                return
            endif

            !Count the number of permutation
            do i=1,size(Indices)
                do j=i,size(Indices)
                    !Find the index that should be at i position
                    !Careful to exclude the case index at i already correct
                    if((j.ne.Indices(j)).and.(Indices(j).eq.i))then
                        !Swap it with the wrong index and count as 1 swap
                        !Repeat until all indices are in place
                        Indices(j)=Indices(i)
                        swapCount=swapCount+1
                        exit
                    endif
                enddo
            enddo

            !Easy converting even swap count to 1 and odd swap count to -1
            e=modulo(swapCount,2)*(-2)+1

        end function LeviCivita

        function modNoZero(i,j)result(m)

            integer :: i,j,m

            m=mod(i,j)
            
            if(m.eq.0)then
                m=j
            endif

        end function modNoZero

end module Functions
