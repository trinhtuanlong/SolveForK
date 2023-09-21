program CubicSolveForK

    use Functions
    use MatrixUtilities
    use Solvers

    complex,parameter :: cEpsilon=(.01000000000000000000000000,0),cTe=(0.0333333333333333333333333333333333,0)

    complex,dimension(3) :: cK

    !complex :: a,b,c,d

    !complex :: p,q,Delta

    ! Root of unity
    !complex,parameter :: epsilon_1=cmplx((-1./2.),(sqrt(3.)/2.)),epsilon_2=cmplx((-1./2.),(-sqrt(3.)/2.))

    ! Eigenvalues
    complex,dimension(3) :: cLambda

    complex,dimension(3,3) :: LHS,V,A

    complex,dimension(3) :: RHS

    ! complex,dimension(3,9) :: LUP

    ! print *, "Starting program"
!
    ! cLambda=(/(.7,0.),(.56,0.),(.28,0.)/)
    ! cLambda=(/(.276,0.),(.276,0.),(.6,0.)/)
    ! cLambda=(/(.47,0.),(.5,0.),(.47,0.)/)
    ! cLambda=(/(.230000000000000000000000,0.),(.2200000000000000000000,0.),(.2100000000000000000000,0.)/)
    cLambda=(/(.99000000000000000000000,0.),(.960000000000000000000,0.),(.93000000000000000000000,0.)/)
    ! cLambda=(/(.70000000000000000000000,0.),(.688000000000000000000,0.),(.6599990000000000000000,0.)/)
    ! cLambda=(/(.670000000000000000000000,0.),(.62000000000000000000,0.),(.5880000000000000000,0.)/)
    ! cLambda=(/(-.04,0.),(.23,0.),(-.04,0.)/)
    ! cLambda=(/(.5396,0.),(.999,0.),(.9277,0.),(.7489,0.)/)
!
    V=IdentityMatrix(3)
    V(size(V,1),1)=-1
    ! do i=1,size(V,1)
    !     print *, V(i,:)
    ! enddo
!
    ! A(1,:)=(/(1.,0.),cTe,(0.,0.)/)
    ! A(2,:)=(/(0.,0.),(1.,0.),(0.,0.)/)
    ! A(3,:)=(/(0.,0.),(0.,0.),-cEpsilon/)
!
    ! A(1,:)=(/(1.,0.),(.02,0.),(0.,0.),(0.,0.)/)
    ! A(2,:)=(/(-.9083,0.),(.6765,0.),(.9083,0.),(0.,0.)/)
    ! A(3,:)=(/(0.,0.),(0.,0.),(1.,0.),(0.,0.)/)
    ! A(4,:)=(/(0.,0.),(0.,0.),(0.,0.),-cEpsilon/)
    ! print *, cLambda
!
    ! LHS(1,:)=(/(5.,0.),(6.,0.),(4.,0.)/)
    ! LHS(2,:)=(/(6.,0.),(4.,0.),(5.,0.)/)
    ! LHS(3,:)=(/(4.,0.),(5.,0.),(6.,0.)/)
!
    ! testDeterminant: block
        ! complex :: detLHS
        ! detLHS=Determinant(LHS)
        ! print *, detLHS
    ! end block testDeterminant

    ! getRHS: block
        ! complex,dimension(1:size(V,1),1:size(A,2)) :: VA
        ! VA=matmul(V,A)
        ! do i=1,size(VA,1)
            ! print *, VA(i,:)
        ! enddo
        ! do i=1,size(RHS,1)
            ! RHS(i)=-CharacteristicPolynomial(VA,cLambda(i))
        ! enddo
        ! print *, RHS
    ! end block getRHS

    ! getLHS: block
        ! complex,dimension(4,4) :: modVA
        ! do i=1,size(cLambda)
            ! modVA=cLambda(i)*IdentityMatrix(4)-matmul(V,A)
            ! print *, (modVA(j,:),NEW_LINE('1'),j=1,size(modVA,2))
            ! modVA(:,1)=modVA(:,1)-modVA(:,size(modVA,2))
            ! modVA(:,size(modVA,2))=0.
            ! print *, (modVA(j,:),NEW_LINE('1'),j=1,size(modVA,2))
            ! do j=1,size(LHS,2)
                ! LHS(i,j)=MatrixMinor(modVA,j,size(modVA,2))
            ! enddo
        ! enddo
        ! do i=1,size(LHS,1)
            ! print *, LHS(i,:)
        ! enddo
    ! end block getLHS
!
    do i=1,size(LHS,1)
        LHS(i,1)=cLambda(i)**2-(2-cEpsilon)*cLambda(i)+(1-cEpsilon)
        LHS(i,2)=cEpsilon*cTe
        LHS(i,3)=cLambda(i)**2-2*cLambda(i)+1
        print *, LHS(i,:)
        RHS(i)=-cLambda(i)**3+(2-cEpsilon)*cLambda(i)**2-(1-2*cEpsilon)*cLambda(i)-cEpsilon
        print *, RHS(i)
    enddo

    ! MaximizeDelta: block
        ! complex :: Delta=(0.,0.),tempDelta
        ! complex,dimension(3) :: cKtemp=(/0.,0.,0./)
        ! complex,dimension(3) :: cLambdatemp
        ! do i1=0,100
            ! do i2=0,100
                ! do i3=0,100
                    ! print *, i1,i2,i3
                    ! cLambdatemp=(/cmplx(.5-(real(i1)*.01),0.),&
                                  ! cmplx(.5-(real(i2)*.01),0.),&
                                  ! cmplx(.5-(real(i3)*.01),0.)/)
                    ! do i=1,size(LHS,1)
                        ! LHS(i,1)=cLambdatemp(i)**2-(2-cEpsilon)*cLambdatemp(i)+(1-cEpsilon)
                        ! LHS(i,2)=cEpsilon*cTe
                        ! LHS(i,3)=cLambdatemp(i)**2-2*cLambdatemp(i)+1
                        ! print *, LHS(i,:)
                        ! RHS(i)=-cLambdatemp(i)**3+(2-cEpsilon)*cLambdatemp(i)**2-(1-2*cEpsilon)*cLambdatemp(i)-cEpsilon
                        ! print *, RHS(i)
                    ! enddo
                    ! SolveLinearBlock: block
                        ! logical :: Solvable
                        ! call SolveLinear(LHS,RHS,cKtemp,Solvable)
                        ! print *, Solvable
                        ! if(Solvable)then
                            ! tempDelta=(6-2*(cLambdatemp(1)+cLambdatemp(2)+cLambdatemp(3)))/&
                                       ! sqrt(cKtemp(1)*conjg(cKtemp(1))+&
                                       ! cKtemp(2)*conjg(cKtemp(2))+&
                                       ! cKtemp(3)*conjg(cKtemp(3)))
                            ! print *, tempDelta
                            ! if(real(tempDelta).gt.real(Delta))then
                                ! Delta=tempDelta
                                ! print *, Delta
                                ! cK=cKtemp
                                ! cLambda=cLambdatemp
                            ! endif
                        ! else
                            ! cycle
                        ! endif
                    ! end block SolveLinearBlock
                ! enddo
            ! enddo
        ! enddo
    ! print *, Delta
    ! end block MaximizeDelta

    testLambda: block
        logical :: Solvable
        call SolveLinear(LHS,RHS,cK,Solvable,argVerbose=.true.)
    end block testLambda
!
    ! print *, cLambda
    print *, cK
    print *, matmul(LHS,cK)-RHS
    ! print *, cK(1)+cK(3)+cEpsilon-2
!
end program CubicSolveForK
