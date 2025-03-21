!
!   Module to estimate drag force on a particle in Stokes flow with respect to neighbouring particles
!
MODULE mNeuralNets

    use mCommon
    
    IMPLICIT NONE

!
!    Parameters
!
! ----------------------------------------------------------------------
!
!    integer, parameter :: rk = selected_real_kind(8)
    integer, parameter :: nnSigmoid = 101
    integer, parameter :: nnRELU = 102
    integer, parameter :: nnIdentity = 103
    integer, parameter :: nnSoftPlus = 104    
    integer, parameter :: nnTanH = 105    
!
!    Types
!
! ----------------------------------------------------------------------
    type Matrix
!
!   Defines the matrix type
!    
        integer nrow,ncol 
        real(rk), pointer :: val(:,:)  ! values of matrix elements
    end type
! ----------------------------------------------------------------------
      type Vector
!
!   Defines the vector type
!    
        integer nrow
        real(rk), pointer :: val(:)
    end type

! ----------------------------------------------------------------------
!
    TYPE NeuralNet
        character(255) name
        integer nLayers ! total number of neruons (including input and output)
        integer nHiddenLayers
        integer nSteps
        integer inpLength,outLength  ! number of rows in first and last layer
        integer, pointer :: AFtype(:)  ! type of activation function        
        integer, pointer :: layerLengths(:)  ! Lengths of layers
        type(Matrix), pointer :: weight(:)  ! weight matrices (nLayer-1)
        type(Vector), pointer :: neuron(:)  ! neurons (nLayer) 
        type(Vector), pointer :: bias(:)
    
    END TYPE

!
! ----------------------------------------------------------------------
!

!
!   Subroutines
!
    CONTAINS
!
! ----------------------------------------------------------------------
!
    Type(NeuralNet) function nnCreate(fileName)
        Type(NeuralNet) nn
        character(*) :: fileName
        character(999) OneLine
        character(255), ALLOCATABLE :: AFtypeString(:)
        integer :: unit, i,r,c

        open(newunit=unit, file=trim(fileName), status='old', action='read', iostat=i)
        if (i /= 0) then
            print *, "Error opening file: ", trim(fileName)
            stop
        end if

        ! Read NN name        
        call rOneTL(unit,OneLine)
        nn%name = TRIM(OneLine)

        ! Read the number of layers
        call rOneTL(unit,OneLine)        
        read(OneLine,*) nn%nLayers
        nn%nSteps = nn%nLayers - 1
        nn%nHiddenLayers = nn%nLayers - 2

        ! Read layer lengths
        allocate (nn%layerLengths(nn%nLayers))
        call rOneTL(unit,OneLine)  
        READ(OneLine,*) (nn%layerLengths(i),i=1,nn%nLayers)

        nn%inpLength = nn%layerLengths(1)
        nn%outLength = nn%layerLengths(nn%nLayers)

        ! Allocate neurons
        allocate(nn%neuron(nn%nLayers))        
        do i = 1, nn%nLayers
            nn%neuron(i)%nrow = nn%layerLengths(i)
            allocate(nn%neuron(i)%val(nn%neuron(i)%nrow))
        end do

        ! Allocate weights
        allocate(nn%weight(nn%nLayers-1))
        do i = 1, nn%nLayers - 1
            nn%weight(i)%nrow = nn%layerLengths(i + 1)
            nn%weight(i)%ncol = nn%layerLengths(i)
            allocate(nn%weight(i)%val(nn%weight(i)%nrow,nn%weight(i)%ncol))
        end do        

        ! Allocate Biases
        allocate(nn%bias(nn%nLayers-1))
        do i = 1, nn%nLayers - 1
            nn%bias(i)%nrow = nn%layerLengths(i + 1)
            allocate(nn%bias(i)%val(nn%bias(i)%nrow))
        end do

        ! Read activation functions
        allocate(nn%AFtype(nn%nLayers-1))
        allocate(AFtypeString(nn%nLayers-1))    
        call rOneTL(unit,OneLine) 
        READ(OneLine,*) (AFtypeString(i),i=1,nn%nLayers-1)

        do i  = 1, nn%nLayers-1
            select case (StrLowCase(TRIM(AFtypeString(i))))
            case ("sigmoid")
                nn%AFtype(i)=nnSigmoid
            case ("identity")
                nn%AFtype(i)=nnIdentity
            case ("relu")
                nn%AFtype(i)=nnRELU    
            case ("tanh")
                nn%AFtype(i)=nnTanH    
            case ("softplus")
                nn%AFtype(i)=nnSoftPlus                                          
            case default
              Print *,"input file ERROR in NeuralNetSetup"
            end select
        end do
        DEALLOCATE(AFtypeString)
!
        ! Read the weight matrices and biases
        do i = 1, nn%nLayers-1
            do r = 1, nn%weight(i)%nrow
                call rOneTL(unit,OneLine) 
                READ(OneLine,*) (nn%weight(i)%val(r, c), c=1,nn%weight(i)%ncol)    
            end do
            do r = 1, nn%bias(i)%nrow
                call rOneTL(unit,OneLine) 
                READ(OneLine,*) nn%bias(i)%val(r)
            end do
        end do

        close(unit)
        nnCreate = nn

    end function

!
! -------------------------------------------------------------------------------------------------
!
    subroutine nnSetNeurons(nn,x)

        Type(NeuralNet) nn
        real(rk) x(nn%inpLength)
        integer i,j
        ! copy x to first neuron
        do i = 1, nn%inpLength
          nn%neuron(1)%val(i) = x(i)
        end do
        do i = 1, nn%nSteps
          nn%neuron(i+1)%val = nn%bias(i)%val + MATMUL(nn%weight(i)%val,nn%neuron(i)%val)
          ! activate
          do j = 1, nn%neuron(i+1)%nrow
            nn%neuron(i+1)%val(j) = nnActivationFunction(nn%neuron(i+1)%val(j),nn%AFtype(i))
          end do
        end do
      
    end subroutine
      
      
!
! -------------------------------------------------------------------------------------------------
!
    subroutine nnEval(nn,x,y)
    
      Type(NeuralNet) nn
      real(rk) x(nn%inpLength)
      real(rk) y(nn%outLength)
      integer i
    
      ! evaluate neurons
      call nnSetNeurons(nn,x)
      ! copy last neuron to y
      do i = 1, nn%outLength
        y(i) = nn%neuron(nn%nLayers)%val(i)
      end do
    
    end subroutine
      

!
! ----------------------------------------------------------------------
!
    real(rk) function nnActivationFunction(x,AFtype)

      real(rk) x,y
      integer AFtype

      select case (AFtype)

        case (nnSigmoid) 
          y = 1.0_rk / (1.0_rk + EXP(-x))

        case (nnRELU)
          if (x<0) then
            y = 0.0_rk
          else
            y = x
          end if

        case (nnIdentity) 
          y = x

        case (nnSoftPlus)
          y = log(1.0_rk + exp(x))

        case (nnTanH)
          y = tanh(x)

        case default
          Print *,"ERROR!"

      end select
      nnActivationFunction = y

    end function



!______________________________________________________________________C
    SUBROUTINE rOneTL(lun,OneLine)
!   _    ___ _    _
!   Read One Text Line
!
!   Returns the first nonempty text line in file LUN, which does not
!   include the # character. If end of file is encoutered, it returns EOF
      CHARACTER*(*) OneLine
      INTEGER lun,i

10    READ(lun,'(A)',END=20) OneLine

!     Check if line is empty
      IF (len_trim(OneLine).EQ.0) GOTO 10

!     Check if line contains # character
      DO i=1,len_trim(OneLine)
        IF (OneLine(i:i).EQ.'#') GOTO 10
      ENDDO

      RETURN

20      OneLine='EOF'
    END
        

    FUNCTION StrLowCase( Input_String ) RESULT( Output_String ) 

        CHARACTER(*), PARAMETER :: LOWER_CASE = 'abcdefghijklmnopqrstuvwxyz' 
        CHARACTER(*), PARAMETER :: UPPER_CASE = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' 

        CHARACTER(*), INTENT(IN)     :: Input_String 
        CHARACTER(LEN(Input_String)) :: Output_String 
        INTEGER :: i, n 

        ! Copy input string 
        Output_String = Input_String 

        ! Convert case character by character 
        DO i = 1, LEN(Output_String) 
          n = INDEX(UPPER_CASE, Output_String(i:i)) 
          IF ( n /= 0 ) Output_String(i:i) = LOWER_CASE(n:n) 
        END DO 
      END FUNCTION StrLowCase     
        

!
! ----------------------------------------------------------------------
!

END MODULE mNeuralNets