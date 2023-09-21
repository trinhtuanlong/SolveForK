module MatrixUtilities

    use Functions

    public Determinant

    interface Determinant
        module procedure DeterminantScalar
        module procedure DeterminantMatrix
    end interface Determinant

    public MatrixMinor
    
    interface MatrixMinor
        module procedure MatrixFirstMinor
    end interface MatrixMinor

    contains

        function IdentityMatrix(pSize)result(I)
            
            integer :: pSize

            integer,dimension(pSize,pSize) :: I

            do j=1,size(I,1)
                do k=1,size(I,2)
                    I(j,k)=KroneckerDelta(j,k)
                enddo
            enddo

        end function IdentityMatrix


        function LUPDecomposition(A)result(LUP)

            complex,dimension(1:,1:) :: A

            complex,dimension(1:size(A,1),1:size(A,2)) :: L

            complex,dimension(1:size(A,1),1:(2*size(A,2))) :: UP

            complex,dimension(1:size(A,1),1:(3*size(A,2))) :: LUP

            logical :: Solvable=.true.

            ! print *, A

            !Check square matrix
            if(size(A,1).ne.size(A,2))then
                print *, "Module MatrixUtilities function LUPDecomposition: Not a square matrix!"
                return
            endif

            !Populate the left half of the UP matrix with the input matrix and the right half with the pivot matrix
            !Populate the L matrix with the identity matrix
            do i=1,size(A,1)
                do j=1,size(A,2)
                    UP(i,j)=A(i,j)
                    UP(i,j+size(A,2))=KroneckerDelta(i,j)
                    L(i,j)=KroneckerDelta(i,j)
                enddo
            enddo

            !Main operations
            do i=1,size(UP,1)
                !Check if the current diagonal element is 0
                if(UP(i,i).eq.(0.,0.))then
                    !Assume unsolvable until a non-zero element is found
                    Solvable=.false.
                    !Check for non-zero element in the same column
                    do j=i,size(UP,1)
                        if(UP(j,i).ne.(0.,0.))then
                            !Set solvable if found
                            Solvable=.true.
                            !Swap the rows
                            swappingRows: block
                                complex,dimension(size(UP,2)) :: tempUP
                                tempUP=UP(i,:)
                                UP(i,:)=UP(j,:)
                                UP(j,:)=tempUP
                            end block swappingRows
                        endif
                    enddo
                    !If no non-zero elements found then return not a square matrix
                    if(.not.Solvable)then
                        print *, "Module MatrixUtilities function LUPDecomposition: Not a square matrix!"
                        return
                    endif
                endif
                !Eliminate columns
                do j=i+1,size(UP,1)
                    !Save the coefficient into the L matrix
                    L(j,i)=UP(j,i)/UP(i,i)
                    !Subtract from each row the current row times the coefficient
                    do k=1,size(A,2)
                        UP(j,k)=UP(j,k)-UP(i,k)*L(j,i)
                        ! do mn=1,size(UP,1)
                            ! print *, UP(mn,:)
                        ! enddo
                    enddo
                enddo
            enddo

            !Save the result into the LUP matrix
            do i=1,size(LUP,1)
                do j=1,size(L,2)
                    LUP(i,j)=L(i,j)
                enddo
                do j=size(L,2)+1,size(LUP,2)
                    LUP(i,j)=UP(i,j-size(L,2))
                enddo
            enddo

            ! print *, (LUP(:,i),NEW_LINE('1'),i=1,size(LUP,2))

        end function LUPDecomposition

        function DeterminantScalar(a)result(cd)
            
            complex :: a

            complex :: cd

            cd=a

        end function DeterminantScalar

        function DeterminantMatrix(A)result(cD)

            complex, dimension(1:,1:) :: A

            complex :: cD,cDetU

            integer :: iDetP

            integer,dimension(size(A,1)) :: PIndices

            complex, dimension(1:size(A,1),3*size(A,2)) :: LUP

            cDetU=1.

            !Get the LUP decomposition
            LUP=LUPDecomposition(A)

            !Get the determinant of the U matrix
            do i=1,size(A,1)
                ! print *, LUP(i,i+size(A,2))
                cDetU=cDetU*LUP(i,i+size(A,2))
                ! print *, cDetU
            enddo
            ! print *, cDetU

            !Get the indices of the unity elements in each column of the P matrix
            do i=1,size(A,1)
                do j=1,size(A,2)
                    if(LUP(i,j+2*size(A,2)).ne.(0.,0.))then
                        PIndices(i)=j
                    endif
                enddo
            enddo

            ! print *, PIndices
            !Put into the Levi-Civita symbol for the determinant
            iDetP=LeviCivita(PIndices)

            ! print *, iDetP

            cD=cDetU/iDetP            

        end function DeterminantMatrix

        function CharacteristicPolynomial(A,lambda)result(p)
            
            complex,dimension(1:,1:) :: A

            complex :: lambda

            complex :: p

            if(size(A,1).ne.size(A,2))then
                print *, "Module MatrixUtilities function CharacteristicPolynomial: Not a square matrix!"
                return
            endif

            ! print *, lambda*IdentityMatrix(size(A,1))-A
            p=Determinant(lambda*IdentityMatrix(size(A,1))-A)

        end function CharacteristicPolynomial

        function MatrixFirstMinor(A,i,j)result(m)

            complex,dimension(1:,1:) :: A

            !i is row index and j is column index
            integer :: i,j

            complex :: m

            !Write the matrix with i row and j column removed here
            complex,dimension(1:(size(A,1)-1),1:(size(A,2)-1)) :: modA

            !First check if indices are within range
            if(i.lt.1.or.j.lt.1)then
                print *, "Module MatrixUtilities function MatrixFirstMinor in MatrixMinor: Indices must be greater than 0"
                return
            elseif(i.gt.size(A,1))then
                print *, "Module MatrixUtilities function MatrixFirstMinor in MatrixMinor: Row index out of row range"
                return
            elseif(j.gt.size(A,2))then
                print *, "Module MatrixUtilities function MatrixFirstMinor in MatrixMinor: Column index out of column range"
                return
            endif

            ! print *,(A(k,:),NEW_LINE('1'),k=1,size(A,1))

            !Start populating the minor matrix
            !First populate the upper left corner
            modA(1:i-1,1:j-1)=A(1:i-1,1:j-1)
            !Only populate the lower left corner if there are rows left
            if(i.lt.size(A,1))then
                modA(i:,1:j-1)=A(i+1:,:j-1)
            endif
            !Only populate the upper right corner if there are columns left
            if(j.lt.size(A,2))then
                modA(:i-1,j:)=A(:i-1,j+1:)
            endif
            !Only populate the lower right corner if there are elements left
            if(i.lt.size(A,1).and.j.lt.size(A,2))then
                modA(i:,j:)=A(i+1:,j+1:)
            endif

            ! print *,(modA(k,:),NEW_LINE('1'),k=1,size(modA,1))

            m=Determinant(modA)

            ! print *, m

        end function MatrixFirstMinor
        
        ! function HouseholderTransformation(A,i)result(B)
        ! end function HouseholderTransformation

end module MatrixUtilities
