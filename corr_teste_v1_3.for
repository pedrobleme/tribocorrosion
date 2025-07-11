!> @brief Módulo para centralizar os parâmetros padrão do Abaqus.
!!
!! Implementação do modelo de corrosão em umeshmotion
!! 2.1. - ok
!! 2.2 - ok

!!v.1.4-> Reescrita de vexternaldb ;SUBROUTINE read_element_connectivity(ok)
!v1.4-v1.5 - Início da implementação de umeshmotion:
!Atualização de vexternaldb
!SUBROUTINE build_node_to_element_connectivity(ok)
!SUBROUTINE UMESHMOTION
!INTEGER FUNCTION FIND_NODE_ID(coords)
!
!--------------------------------------------------------------------------------------------------------------
      MODULE AbaqusParameters
      IMPLICIT NONE

    ! --- Parâmetros de Controle do Job (lOp) ---
      INTEGER, PARAMETER :: J_INT_START_ANALYSIS    = 0
      INTEGER, PARAMETER :: J_INT_START_STEP        = 1
      INTEGER, PARAMETER :: J_INT_SETUP_INCREMENT   = 2
      INTEGER, PARAMETER :: J_INT_START_INCREMENT   = 3
      INTEGER, PARAMETER :: J_INT_END_INCREMENT     = 4
      INTEGER, PARAMETER :: J_INT_END_STEP          = 5
      INTEGER, PARAMETER :: J_INT_END_ANALYSIS      = 6

    ! --- Índices para o i_Array (Informações do Job) ---
      INTEGER, PARAMETER :: I_INT_NTOTALNODES       = 1
      INTEGER, PARAMETER :: I_INT_NTOTALELEMENTS    = 2
      INTEGER, PARAMETER :: I_INT_KSTEP             = 3
      INTEGER, PARAMETER :: I_INT_KINC              = 4
      INTEGER, PARAMETER :: I_INT_ISTATUS           = 5
      INTEGER, PARAMETER :: I_INT_LWRITERESTART     = 6

    ! --- Índices para o r_Array (Informações de Tempo) ---
      INTEGER, PARAMETER :: I_FLT_TOTALTIME         = 1
      INTEGER, PARAMETER :: I_FLT_STEPTIME          = 2
      INTEGER, PARAMETER :: I_FLT_DTIME             = 3

      END MODULE AbaqusParameters
!--------------------------------------------------------------------------------------------------------------
!> @brief Declara variáveis globais, parâmetros e matrizes compartilhadas.
!!
!! Este módulo centraliza todas as variáveis que precisam ser acessadas por diferentes sub-rotinas da simulação.
!--------------------------------------------------------------------------------------------------------------
      MODULE SharedVariables
      IMPLICIT NONE

    ! --- Precisão Numérica ---
      INTEGER, PARAMETER :: dp = KIND(1.0d0)

    ! --- Parâmetros da Simulação ---
      INTEGER             :: n_elems             = 0    ! Número total de elementos na malha.
      INTEGER, PARAMETER  :: max_influence_elems = 500  ! Número máximo de elementos no raio de influência.
      INTEGER, PARAMETER  :: max_face_neighbors = 6

    ! --- Parâmetros Físicos e do Modelo ---
      REAL(KIND=dp), PARAMETER :: lr    = 0.3      ! [mm] Comprimento intrínseco do modelo não-local.
      REAL(KIND=dp), PARAMETER :: nt    = 10950.0  ! [dias] Tempo total de corrosão.
      REAL(KIND=dp), PARAMETER :: sigth = 121      ! [MPa] Tensão de escoamento do material.

    ! --- Índices para a Matriz elem_properties ---
    ! Esta matriz armazena todas as propriedades lidas do arquivo e informações de topologia para cada elemento.
      INTEGER, PARAMETER :: PROP_IDX_ELEM_ID           = 1   ! Coluna 1: Número do elemento (ID)
      INTEGER, PARAMETER :: PROP_IDX_NEIGHBOR_1        = 2   ! Coluna 2: ID do vizinho na face 1
      INTEGER, PARAMETER :: PROP_IDX_NEIGHBOR_2        = 3   ! Coluna 3: ID do vizinho na face 2
      INTEGER, PARAMETER :: PROP_IDX_NEIGHBOR_3        = 4   ! Coluna 4: ID do vizinho na face 3
      INTEGER, PARAMETER :: PROP_IDX_NEIGHBOR_4        = 5   ! Coluna 5: ID do vizinho na face 4
      INTEGER, PARAMETER :: PROP_IDX_NEIGHBOR_5        = 6   ! Coluna 6: ID do vizinho na face 5
      INTEGER, PARAMETER :: PROP_IDX_NEIGHBOR_6        = 7   ! Coluna 7: ID do vizinho na face 6
      INTEGER, PARAMETER :: PROP_IDX_SURFACE_FLAG      = 8   ! Coluna 8: Flag que indica se o elemento é de superfície (1.0) ou não (0.0)
      INTEGER, PARAMETER :: PROP_IDX_PITTING           = 9   ! Coluna 9: Valor de pitting normalizado
      INTEGER, PARAMETER :: PROP_IDX_VOLUME            = 10  ! Coluna 10: Volume do elemento
      INTEGER, PARAMETER :: PROP_IDX_COORD_X           = 11  ! Coluna 11: Coordenada X do centroide
      INTEGER, PARAMETER :: PROP_IDX_COORD_Y           = 12  ! Coluna 12: Coordenada Y do centroide
      INTEGER, PARAMETER :: PROP_IDX_COORD_Z           = 13  ! Coluna 13: Coordenada Z do centroide
      INTEGER, PARAMETER :: PROP_IDX_PITTING_ABSOLUTE  = 14  ! Coluna 14: Valor de pitting "absolu 

    !-----------------------------------------------------------------------
    ! --- Flags e Estados de Controle da Análise ---
    !-----------------------------------------------------------------------
      LOGICAL :: failure_occurred_this_increment = .FALSE.

    ! --- Marca o Início da Corrosão ---
      INTEGER, PARAMETER :: flag_corrosion_started = 10

    ! --- Controle de Deleção de Elemento ---
      INTEGER, PARAMETER :: delete_element_flag = 0      ! Flag para ativar a deleção de elemento via stateNew(k,3).

    ! --- Status de Corrosão do Elemento (element_corrosion_status) ---
      INTEGER, PARAMETER :: flag_corrosion_active   = 0    ! O elemento está sujeito à corrosão.
      INTEGER, PARAMETER :: flag_corrosion_inactive = -10  ! O elemento não sofre mais corrosão.

    ! --- Status de Deleção por Deformação (strain_deletion_status) ---
      INTEGER, PARAMETER :: flag_deletion_by_strain_enabled  = 0    ! A deleção por deformação está ativa.
      INTEGER, PARAMETER :: flag_deletion_by_strain_disabled = -20  ! A deleção por deformação está inativa.

    ! --- Trava de Atualização de Propriedades (property_update_lock) ---
      INTEGER, PARAMETER :: key_prop_update_unlocked = 0    ! As propriedades do elemento podem ser atualizadas.
      INTEGER, PARAMETER :: key_prop_update_locked   = -5   ! As propriedades do elemento estão travadas.

    ! --- Gatilho para Transferência Geral de Propriedades ---
      INTEGER, PARAMETER :: flag_general_property_transfer        = -8
      INTEGER, PARAMETER :: flag_general_property_transfer_locked = 0

    ! --- Flag para Parada Global da Corrosão (global_corrosion_stop_flag) ---
      INTEGER, PARAMETER :: flag_global_stop = -15

    ! --- Flags de Controle de Execução ---
      INTEGER :: prop_update_trigger         ! Gatilho para a rotina transfer_properties.
      INTEGER :: global_corrosion_stop_flag  ! Controla a parada global do processo de corrosão.

    ! --- Propriedades do Material e Simulação ---
      REAL(KIND=dp)                                 :: mCO, timetot, lmCO
      REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE    :: elem_properties

    ! --- Matrizes e Vetores Compartilhados ---
      INTEGER,       DIMENSION(:),   ALLOCATABLE    :: element_corrosion_status
      INTEGER,       DIMENSION(:),   ALLOCATABLE    :: strain_deletion_status
      INTEGER,       DIMENSION(:),   ALLOCATABLE    :: property_update_lock
      INTEGER,       DIMENSION(:),   ALLOCATABLE    :: counter
      INTEGER,       DIMENSION(:,:), ALLOCATABLE    :: influ
      REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE    :: dist
      REAL(KIND=dp), DIMENSION(:),   ALLOCATABLE    :: temp_pitting_write
      REAL(KIND=dp), DIMENSION(:),   ALLOCATABLE    :: temp_surface_flag_write
      REAL(KIND=dp), DIMENSION(:),   ALLOCATABLE    :: nonlocal_property
      REAL(KIND=dp), DIMENSION(:),   ALLOCATABLE    :: temp_nonlocal_write
      REAL(KIND=dp), DIMENSION(:),   ALLOCATABLE    :: damage_current
      REAL(KIND=dp), DIMENSION(:),   ALLOCATABLE    :: von_mises_stress
      REAL(KIND=dp), DIMENSION(:),   ALLOCATABLE    :: non_local_stress
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE    :: element_recession_velocity
      REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE    :: initial_nodal_normals !(n_total_nodes, 3)

    ! --- Variáveis para Tribocorrosão ---      
      REAL(KIND=dp) :: accumulated_pure_tribo
      REAL(KIND=dp) :: enhanced_tribo_damage_inc
      REAL(KIND=dp) :: combined_damage, total_damage
      REAL(KIND=dp) :: corrosion_depth_inc
      REAL(KIND=dp) :: v_cor
      REAL(KIND=dp) :: wear_depth_inc
      REAL(KIND=dp) :: v_wear
      REAL(KIND=dp) :: v_total_scalar
      REAL(KIND=dp) :: element_recession_vector(3)
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE    :: tribological_damage_accumulated

      INTEGER :: n_nodes ! Número total de nós no modelo !FIX
      REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: nodal_normals
      !REAL(KIND=dp), DIMENSION(:), ALLOCATABLE   :: relative_slip_increment
      !REAL(KIND=dp), DIMENSION(:), ALLOCATABLE   :: contact_pressure

      !!-------- v1.3->v1.4

      INTEGER, PARAMETER :: max_nodes_per_elem = 8  ! Para elementos C3D8
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: element_nodes  ! [n_elems, max_nodes_per_elem]
      INTEGER, DIMENSION(:), ALLOCATABLE :: num_nodes_per_element ! [n_elems]

      !----v1.4->v1.5----

      ! --- Variáveis para UMESHMOTION ---
      INTEGER, DIMENSION(:), ALLOCATABLE :: node_elements_count       ! Número de elementos conectados a cada nó
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: node_to_elem_conn       ! Conectividade nó-elemento
      REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: nodal_velocities  ! Velocidades nodais calculadas (para UMESHMOTION)
      REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: nodal_coordinates
      
      CONTAINS
      INTEGER FUNCTION FIND_CLOSEST_NODE(coords)
          REAL(KIND=dp), INTENT(IN) :: coords(3)
          REAL(KIND=dp) :: min_dist, distance, node_coords(3)
          INTEGER :: i
          
          min_dist = HUGE(1.0_dp)
          FIND_CLOSEST_NODE = 0
          
          IF (.NOT. ALLOCATED(nodal_coordinates)) THEN
              WRITE(*,*) 'ERRO: nodal_coordinates não alocado em FIND_CLOSEST_NODE'
              RETURN
          END IF
          
          DO i = 1, SIZE(nodal_coordinates,1)
              node_coords = nodal_coordinates(i,:)
              distance = SQRT(SUM((node_coords - coords)**2))
              IF (distance < min_dist) THEN
                  min_dist = distance
                  FIND_CLOSEST_NODE = i
              END IF
              
              IF (min_dist < 1.0e-6_dp) EXIT
          END DO
          
          IF (min_dist > 1.0_dp) THEN
              WRITE(*,*) 'AVISO: Nó distante encontrado - distance:', min_dist
          END IF
      END FUNCTION FIND_CLOSEST_NODE

      END MODULE SharedVariables
!-------------------------------------------------------------------------------------------------------------------------------------
!> @brief Ponto de entrada principal da VUMAT chamado pelo Abaqus.
!!   
! Esta sub-rotina serve como um "wrapper" ou "invólucro". Sua função é desempacotar os argumentos DO vetor 'jblock'
! (como número DO elemento e ponto de integração) e passá-los de forma organizada para a sub-rotina de trabalho principal, vumatXtrArg.
!-------------------------------------------------------------------------------------------------------------------------------------
      SUBROUTINE vumat (
! Read only -
     *     jblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     *     stepTime, totalTime, dt, cmname, coordMp, charLength,
     *     props, density, strainInc, relSpinInc,
     *     tempOld, stretchOld, defgradOld, fieldOld,
     *     stressOld, stateOld, enerInternOld, enerInelasOld,
     *     tempNew, stretchNew, defgradNew, fieldNew,
! Write only -
     *     stressNew, stateNew, enerInternNew, enerInelasNew )
!
      INCLUDE 'vaba_param.inc'
!     
!      USE subs
      DIMENSION jblock(*), props(nprops),density(*), coordMp(*),
     1     charLength(*), strainInc(*),
     2     relSpinInc(*), tempOld(*),
     3     stretchOld(*),
     4     defgradOld(*),
     5     fieldOld(*), stressOld(*),
     6     stateOld(*), enerInternOld(*),
     7     enerInelasOld(*), tempNew(*),
     8     stretchNew(*),
     9     defgradNew(*),
     1     fieldNew(*),
     2     stressNew(*), stateNew(*),
     3     enerInternNew(*), enerInelasNew(*)

      character*80 cmname
      PARAMETER (i_umt_nblock = 1,
     *     i_umt_npt    = 2,
     *     i_umt_layer  = 3,
     *     i_umt_kspt   = 4,
     *     i_umt_noel   = 5 )

      CALL  vumatXtrArg ( jblock(i_umt_nblock),
     *     ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     *     stepTime, totalTime, dt, cmname, coordMp, charLength,
     *     props, density, strainInc, relSpinInc,
     *     tempOld, stretchOld, defgradOld, fieldOld,
     *     stressOld, stateOld, enerInternOld, enerInelasOld,
     *     tempNew, stretchNew, defgradNew, fieldNew,
     *     stressNew, stateNew, enerInternNew, enerInelasNew,
     *     jblock(i_umt_noel), jblock(i_umt_npt),
     *     jblock(i_umt_layer), jblock(i_umt_kspt))

      RETURN
      END SUBROUTINE vumat
!-------------------------------------------------------------------------------------------------------------------------------------
!> @brief Implementação principal da lógica do modelo constitutivo.
!!
!! Esta sub-rotina contém todo o cálculo do modelo de dano (corrosão,!! desgaste, tribocorrosão)
!! e plasticidade. Ela recebe uma lista de!! argumentos mais clara da rotina 'vumat' para facilitar a leitura e a manutenção do código.
!-------------------------------------------------------------------------------------------------------------------------------------  
      SUBROUTINE vumatXtrArg (
  ! Argumentos de entrada (Read only)
     *     nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     *     stepTime, totalTime, dt, cmname, coordMp, charLength,
     *     props, density, strainInc, relS  pinInc,
     *     tempOld, stretchOld, defgradOld, fieldOld,
     *     stressOld, stateOld, enerInternOld, enerInelasOld,
     *     tempNew, stretchNew, defgradNew, fieldNew,
  ! Argumentos de saída (Write only)
     *     stressNew, stateNew, enerInternNew, enerInelasNew,
  ! Argumentos extras (Read only)
     *     nElement, nMatPoint, nLayer, nSecPoint )
     
      USE SharedVariables
      INCLUDE 'vaba_param.inc'
  !-----------------------------------------------------------------------
  ! --- 1. DECLARAÇÃO DE VARIÁVEIS E ARGUMENTOS ---
  !-----------------------------------------------------------------------
    ! --- Argumentos da Sub-rotina ---
      DIMENSION props(nprops), density(nblock), coordMp(nblock,*),
     1     charLength(nblock), strainInc(nblock,ndir+nshr),
     2     relSpinInc(nblock,nshr), tempOld(nblock),
     3     stretchOld(nblock,ndir+nshr),
     4     defgradOld(nblock,ndir+nshr+nshr),
     5     fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     6     stateOld(nblock,nstatev), enerInternOld(nblock),
     7     enerInelasOld(nblock), tempNew(nblock),
     8     stretchNew(nblock,ndir+nshr),
     9     defgradNew(nblock,ndir+nshr+nshr),
     1     fieldNew(nblock,nfieldv),
     2     stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     3     enerInternNew(nblock), enerInelasNew(nblock)
    ! Documentation of extra arguments:
    !  nElement: Array of internal element numbers
    !  nMatPoint: Integration point number
    !  nLayer   : Layer number for composite shells and layered solids
    !  nSecPoint: Section point number within the current layer
      DIMENSION nElement(nblock), nLayer(nblock),nMatPoint(nblock)

          character*80 cmname
      ! --- Variáveis de Cálculo Locais ---
        INTEGER :: k, nvalue, OK

      ! Variáveis de trabalho para tensões e deformações
        REAL(KIND=dp) :: s11, s22, s33, s12, s13, s23, trace
        REAL(KIND=dp) :: smean, sigdif, vmises, factor
        REAL(KIND=dp) :: peeqOld_k, deqps, eqps, plasticWorkInc
        REAL(KIND=dp) :: yieldOld, yieldNew, hard
        REAL(kind=dp) :: e, xnu, twomu, alamda, thremu
        REAL(kind=dp) ::  masslimit, deltat !Verificar a necessidade dessas duas linhas

      ! Variáveis de trabalho para energia
        REAL(KIND=dp) :: stressPower

      ! Variáveis de trabalho para corrosão e dano
        REAL(KIND=dp) :: tnt, ku, corrosion_damage
        REAL(KIND=dp) :: dCOld_k, dSCld_k
        REAL(KIND=dp) :: sigmaJK, epsilonJK, epsilonf
        REAL(KIND=dp) :: initial_volume
      
      ! Variáveis de trabalho para tribologia
        REAL(KIND=dp) :: calculated_slip_inc, contact_area
        REAL(KIND=dp) :: force_of_friction, energy_dissipated_increment
        REAL(KIND=dp) :: wear_volume_increment, tribological_damage
        REAL(KIND=dp) :: tribological_damage_inc

      ! Arrays Locais
        REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: damage_new
        REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: relative_slip_increment
        REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: contact_pressure
        REAL(KIND=dp), DIMENSION(n_elems)        :: elem_length
        REAL(KIND=dp), DIMENSION(n_elems)        :: corroded_volume 

  !-----------------------------------------------------------------------
  ! --- 2. ATRIBUIÇÃO DE PARÂMETROS FÍSICOS E DE MODELAGEM ---
  !-----------------------------------------------------------------------
      ! --- Constantes Matemáticas e Numéricas ---
      REAL(KIND=dp), PARAMETER :: zero         = 0.0_dp
      REAL(KIND=dp), PARAMETER :: one          = 1.0_dp
      REAL(KIND=dp), PARAMETER :: two          = 2.0_dp
      REAL(KIND=dp), PARAMETER :: half         = 0.5_dp
      REAL(KIND=dp), PARAMETER :: third        = 1.0_dp / 3.0_dp
      REAL(KIND=dp), PARAMETER :: one_point_five = 1.5_dp

      ! --- Parâmetros do Modelo de Dano Geral ---
      REAL(KIND=dp), PARAMETER :: dmax = 0.999_dp         ! Valor máximo do dano antes da deleção do elemento [adimensional]
      REAL(KIND=dp), PARAMETER :: beta = 0.8_dp           ! Fator de propagação de dano para vizinhos [adimensional]

      ! --- Parâmetros do Modelo de Corrosão ---
      REAL(KIND=dp), PARAMETER :: corrosion_rate       = 0.1417_dp ! Taxa de corrosão [mm/ano]
      REAL(KIND=dp), PARAMETER :: corrosion_delay_time = 0.0_dp    ! Tempo de atraso para o início da corrosão [dias]

      ! --- Parâmetros do Modelo de Dano (Johnson-Cook) ---
      REAL(KIND=dp), PARAMETER :: D1 = 0.05_dp  ! Coeficiente D1 do modelo de falha de Johnson-Cook [adimensional]
      REAL(KIND=dp), PARAMETER :: D2 = 3.44_dp  ! Coeficiente D2 do modelo de falha de Johnson-Cook [adimensional]
      REAL(KIND=dp), PARAMETER :: D3 = -2.12_dp ! Coeficiente D3 do modelo de falha de Johnson-Cook [adimensional]
      REAL(KIND=dp), PARAMETER :: D4 = 0.002_dp ! Coeficiente D4 do modelo de falha de Johnson-Cook [adimensional]

      ! --- Parâmetros do Modelo de Desgaste (Tribologia) ---
      REAL(KIND=dp), PARAMETER :: slips_per_day      = 0.10_dp    ! [ciclos/dia] Frequência de deslizamento
      REAL(KIND=dp), PARAMETER :: slip_amplitude     = 1.0_dp     ! [mm] Amplitude do deslizamento
      REAL(KIND=dp), PARAMETER :: test_pressure      = 50.0_dp    ! [MPa] Pressão de contato de referência
      REAL(KIND=dp), PARAMETER :: friction_coefficient = 0.4_dp   ! Coeficiente de atrito [adimensional]
      REAL(KIND=dp), PARAMETER :: wear_coefficient   = 1.0e-5_dp  ! Coef. de desgaste [mm^3/(N*mm)]

      ! --- Parâmetros de Interação (Tribocorrosão) ---
      REAL(KIND=dp), PARAMETER :: coupling_factor_alpha = 0.5_dp ! Fator de sinergia [adimensional]

      ! --- Propriedades Físicas e de Controle ---
      REAL(KIND=dp), PARAMETER :: mass        = 7.86e-6_dp ! [ton/mm^3] Densidade do material
      REAL(KIND=dp), PARAMETER :: corr        = 10.0_dp    ! Fator para critério de massa corroída 

!-------------------------------------------------------------------------------------------------------------------------------------
  ! --- 3. DECLARAÇÃO DE ÍNDICES PARA VARIÁVEIS DE ESTADO (SDVs) ---
  ! SDVs (Solution-Dependent State Variables) são usadas para armazenar o estado do material em cada ponto de integração.
!-------------------------------------------------------------------------------------------------------------------------------------
      INTEGER, PARAMETER :: SDV_PEEQ                 = 1  ! Deformação plástica equivalente
      INTEGER, PARAMETER :: SDV_TOTAL_DAMAGE         = 2  ! Dano total acumulado
      INTEGER, PARAMETER :: SDV_DELETE_FLAG          = 3  ! Flag para deleção do elemento
      INTEGER, PARAMETER :: SDV_PITTING_FACTOR       = 4  ! Fator de pitting
      INTEGER, PARAMETER :: SDV_CORROSION_STATUS     = 5  ! Status da corrosão (ativo/inativo)
      INTEGER, PARAMETER :: SDV_NONLOCAL_PROPERTY    = 6  ! Propriedade não-local
      INTEGER, PARAMETER :: SDV_STRESS_S11           = 7  ! Componente de tensão S11
      INTEGER, PARAMETER :: SDV_STRESS_S22           = 8  ! Componente de tensão S22
      INTEGER, PARAMETER :: SDV_STRESS_S33           = 9  ! Componente de tensão S33
      INTEGER, PARAMETER :: SDV_STRESS_S12           = 10 ! Componente de tensão S12
      INTEGER, PARAMETER :: SDV_STRESS_S13           = 11 ! Componente de tensão S13
      INTEGER, PARAMETER :: SDV_STRESS_S23           = 12 ! Componente de tensão S23
      INTEGER, PARAMETER :: SDV_VON_MISES            = 13 ! Tensão de Von Mises
      INTEGER, PARAMETER :: SDV_DCOld                = 14 ! (Não utilizado, mantido por compatibilidade)
      INTEGER, PARAMETER :: SDV_DSCLd                = 15 ! (Não utilizado, mantido por compatibilidade)
      INTEGER, PARAMETER :: SDV_NONLOCAL_STRESS      = 16 ! Tensão não-local
      INTEGER, PARAMETER :: SDV_SURFACE_FLAG         = 17 ! Flag de superfície
      INTEGER, PARAMETER :: SDV_TOTAL_CORR_DAMAGE    = 18 ! Dano total por corrosão
      INTEGER, PARAMETER :: SDV_CURRENT_TOTAL_DAMAGE = 19 ! Dano total corrente
      INTEGER, PARAMETER :: SDV_JC_DAMAGE_INC        = 20 ! Incremento de dano de Johnson-Cook
      INTEGER, PARAMETER :: SDV_CORROSION_DEPTH      = 21 ! Profundidade da corrosão
      INTEGER, PARAMETER :: SDV_CORROSION_START_TIME = 22 ! Tempo de início da corrosão
      INTEGER, PARAMETER :: SDV_CORROSION_INIT_FLAG  = 24 ! Flag de inicialização da corrosão
      INTEGER, PARAMETER :: SDV_TRIBO_DAMAGE_ACCUM   = 26 ! Dano acumulado por tribologia
      INTEGER, PARAMETER :: SDV_CORRODED_VOLUME      = 30 ! Volume corroído
      
  !-----------------------------------------------------------------------
  ! --- 4. ALOCAÇÃO DE ARRAYS DINÂMICOS ---
  !-----------------------------------------------------------------------
      ALLOCATE (damage_new(n_elems))
      ALLOCATE (relative_slip_increment(n_elems))
      ALLOCATE (contact_pressure(n_elems))

  !-----------------------------------------------------------------------
  ! --- 5. OBTENÇÃO DAS PROPRIEDADES ELÁSTICAS ---
  ! props(1): Módulo de Young, props(2): Coeficiente de Poisson
  ! props(3..) - Dados de escoamento e encruamento
  !-----------------------------------------------------------------------
      e         = props(1)
      xnu       = props(2)
      twomu     = e / (one + xnu)
      alamda    = xnu * twomu / (one - (two * xnu))
      thremu    = one_point_five * twomu
      nvalue    = (nprops / 2) - 1
      masslimit = mass * corr

  !-----------------------------------------------------------------------
  ! --- 6. LÓGICA PRINCIPAL DO MATERIAL ---
  !-----------------------------------------------------------------------
      IF (totalTime .eq. zero) THEN
        ! --- Chute Inicial Elástico (Primeiro Incremento) ---
        ! Essencial para a estabilidade inicial em simulações explícitas.
        DO k = 1, nblock
            trace = strainInc(k,1) + strainInc(k,2) + strainInc(k,3)
            stressNew(k,1) = stressOld(k,1) + twomu * strainInc(k,1) + alamda * trace
            stressNew(k,2) = stressOld(k,2) + twomu * strainInc(k,2) + alamda * trace
            stressNew(k,3) = stressOld(k,3) + twomu * strainInc(k,3) + alamda * trace
            stressNew(k,4) = stressOld(k,4) + twomu * strainInc(k,4)
            IF (nshr > 1) THEN
                stressNew(k,5) = stressOld(k,5) + twomu * strainInc(k,5)
                stressNew(k,6) = stressOld(k,6) + twomu * strainInc(k,6)
            END IF
        END DO
        lmCO = 0.0

      ELSE
        ! --- Início dos Cálculos para Incrementos Posteriores ---
        deltat = nt * dt

        ! --- Laço Principal sobre os Pontos de Integração ---
        DO k = 1, nblock

            ! --- Inicialização de variáveis do laço ---
            corrosion_damage      = 0.0_dp
            tribological_damage   = 0.0_dp
            tribological_damage_inc = 0.0_dp
            deqps                 = 0.0_dp
            s13                   = 0.0_dp
            s23                   = 0.0_dp
            wear_depth_inc        = 0.0_dp


            ! Recupera valores do incremento anterior das variáveis de estado
            peeqOld_k = stateOld(k, SDV_PEEQ)
            dCOld_k   = stateOld(k, SDV_DCOld)
            dSCld_k   = stateOld(k, SDV_DSCLd)
            damage_new(nElement(k)) = stateOld(k, SDV_TOTAL_DAMAGE)
      

  !===============================================================
  ! --- 7. CÁLCULO DE DANO (CORROSÃO, DESGASTE, TRIBOCORROSÃO) ---
  !===============================================================
      
        IF ((totalTime < one) .AND. (element_corrosion_status(nElement(k)) /= flag_corrosion_inactive)) THEN
    !-----------------------------------------------------------
    ! --- 7.1. DANO "PURO" POR CORROSÃO ---
    !-----------------------------------------------------------
          ! tnt: tempo ajustado considerando atraso (corrosion_delay_time)
          tnt = (totalTime-corrosion_delay_time+dt) * nt                 
          ku = 1.0 ! ku: parâmetro cinético de corrosão (ajustável)      
          ! ku = 0.3125 * (tnt**(-0.68))  ! Ex: Valor calibrado experimentalmente
         
        ! Se for o primeiro incremento de corrosão, armazena o tempo inicial
        IF (stateOld(k, SDV_CORROSION_INIT_FLAG) /= flag_corrosion_started) THEN
            stateNew(k, SDV_CORROSION_INIT_FLAG) = flag_corrosion_started
            stateNew(k, SDV_CORROSION_START_TIME) = tnt
        END IF
        
        ! Calcula a profundidade da corrosão                                                   
        stateNew(k, SDV_CORROSION_DEPTH) = ABS(corrosion_rate * (tnt - stateNew(k, SDV_CORROSION_START_TIME)) / 365.0)
     1  * elem_properties(nElement(k),PROP_IDX_SURFACE_FLAG)

        IF (stateNew(k, SDV_CORROSION_DEPTH) < stateOld(k, SDV_CORROSION_DEPTH)) THEN
            stateNew(k, SDV_CORROSION_DEPTH) = stateOld(k, SDV_CORROSION_DEPTH)
        END IF

        ! Calcula o volume corroído e o dano correspondente
        elem_length(nElement(k)) = elem_properties(nElement(k), PROP_IDX_VOLUME) ** (1.0/3.0)             
        corroded_volume(nElement(k)) = ABS(elem_length(nElement(k)) * elem_length(nElement(k))
     1  * stateNew(k, SDV_CORROSION_DEPTH)) * nonlocal_property(nElement(k))
         
         
        initial_volume = elem_properties(nElement(k), PROP_IDX_VOLUME)

        IF (initial_volume > 0.0_dp) THEN
          IF (corroded_volume(nElement(k)) >= initial_volume) THEN
            corrosion_damage = dmax
            stateNew(k, SDV_CORROSION_DEPTH) = 0.0_dp
          ELSE
            corrosion_damage = (ABS(corroded_volume(nElement(k))) * ku) / initial_volume
          END IF
        ELSE
        corrosion_damage = 0.0_dp
        END IF

        !Cálculo da velocidade de avanço da corrosão
        corrosion_depth_inc = stateNew(k,SDV_CORROSION_DEPTH) - stateOld(k,SDV_CORROSION_DEPTH)
        IF (dt > 0.0_dp) THEN
          v_corr = corrosion_depth_inc / dt
        ELSE
          v_corr = 0.0_dp
        END IF


    !-----------------------------------------------------------
    ! --- 7.2. DANO "PURO" POR DESGASTE ---
    !-----------------------------------------------------------
    ! Modelo conceitual de desgaste baseado em energia
        calculated_slip_inc = slips_per_day * slip_amplitude * deltat
        relative_slip_increment(nElement(k)) = calculated_slip_inc
        contact_pressure(nElement(k)) = test_pressure

        IF (relative_slip_increment(nElement(k)) > zero) THEN
        contact_area = elem_properties(nElement(k), PROP_IDX_VOLUME)**(2.0/3.0)
        force_of_friction = contact_pressure(nElement(k)) * friction_coefficient * contact_area
        energy_dissipated_increment = force_of_friction * relative_slip_increment(nElement(k))
        wear_volume_increment = wear_coefficient * energy_dissipated_increment

          ! Cálculo da velocidade de desgast
          IF (contact_area > 0.0_dp) THEN
              wear_depth_inc = wear_volume_increment / contact_area
          END IF

          IF (dt > 0.0_dp) THEN
              v_wear = wear_depth_inc / dt
          ELSE
              v_wear = 0.0_dp
          END IF

          ! Converte volume de desgaste em dano incremental
          IF (elem_properties(nElement(k), PROP_IDX_VOLUME) > zero) THEN
                        tribological_damage_inc = wear_volume_increment / elem_properties(nElement(k), PROP_IDX_VOLUME)            
          ELSE
                tribological_damage = zero
          END IF
        END IF
    !-----------------------------------------------------------
    ! --- 7.3. VELOCIDADE DE RECUO COMBINADA (PARA UMESHMOTION) ---
    !-----------------------------------------------------------
      ! 1. Soma as magnitudes das velocidades de corrosão e desgaste
      v_total_scalar = v_corr + v_wear
      element_recession_velocity(nElement(k)) = v_total_scalar


      ! 2. Define o VETOR de recuo.
      !    Para começar, usamos a simplificação de que o recuo ocorre sempre
      !    na direção -Z global. O ideal seria usar o vetor normal à superfície.
      !element_recession_vector(1) = 0.0_dp  ! Componente X da velocidade
      !element_recession_vector(2) = 0.0_dp  ! Componente Y da velocidade
      !element_recession_vector(3) = -v_total_scalar ! Componente Z da velocidade

      ! 3. Armazena o vetor de velocidade do elemento no array compartilhado.
      !    Este array será lido pela VEXTERNALDB para calcular as velocidades nodais.
      !IF (ALLOCATED(element_recession_velocity)) THEN
      !    element_recession_velocity(nElement(k), :) = element_recession_vector
      !END IF
    !-----------------------------------------------------------
    ! --- 7.4. DANO COMBINADO (TRIBOCORROSÃO) ---
    !-----------------------------------------------------------
        ! Recupera o dano por desgaste acumulado do passo anterior
        accumulated_pure_tribo = stateOld(k, SDV_TRIBO_DAMAGE_ACCUM)

        ! Sinergia: corrosão acelera o novo incremento de desgaste
        enhanced_tribo_damage_inc = tribological_damage_inc * (one + coupling_factor_alpha * corrosion_damage)
                
        ! Acumula o dano por desgaste (já com efeito sinérgico)
        accumulated_pure_tribo = accumulated_pure_tribo + enhanced_tribo_damage_inc

        ! Dano total é a soma da corrosão pura + desgaste acumulado (com sinergia)
        combined_damage = corrosion_damage + accumulated_pure_tribo
        total_damage = MIN(combined_damage, dmax)
        damage_new(nElement(k)) = total_damage
        timetot = totalTime - corrosion_delay_time

  !-----------------------------------------------------------
  ! --- 8. VERIFICAÇÃO DE MASSA E FALHA DO ELEMENTO POR DANO ---
  !-----------------------------------------------------------
    ! Condição de parada por massa corroída
      IF (mCO >= masslimit) THEN
        !write(*,*) mCO, timetot
        !write(*,*)'Fim da Corrosao 100% DO volume atingido'
        CALL block_all(ok)
      END IF

    ! Condição de falha do elemento por dano máximo
        IF (damage_new(nElement(k)) >= dmax) THEN
          stateNew(k, SDV_DELETE_FLAG) = delete_element_flag
          damage_new(nElement(k))      = dmax
          element_corrosion_status(nElement(k)) = flag_corrosion_inactive
          prop_update_trigger = flag_general_property_transfer
          CALL process_failed_element_neighbors(nElement(k), beta, OK)
          failure_occurred_this_increment = .TRUE.
          stateNew(k, SDV_CORROSION_DEPTH) = 0.0
        END IF
      END IF ! Fim do bloco de cálculo de dano

    ! Garante que o dano nunca decresça
      IF (damage_new(nElement(k)) < stateOld(k, SDV_TOTAL_DAMAGE)) THEN
          damage_new(nElement(k)) = stateOld(k, SDV_TOTAL_DAMAGE)
      END IF

  !===============================================================
  ! --- 9. CÁLCULO DE PLASTICIDADE (CRITÉRIO DE VON MISES) ---
  !===============================================================
      CALL vuhard(yieldOld, hard, peeqOld_k, props(3), nvalue)

    ! Tensão de teste elástica
      trace = strainInc(k,1) + strainInc(k,2) + strainInc(k,3)
      s11 = stateOld(k, SDV_STRESS_S11) + twomu * strainInc(k,1) + alamda * trace
      s22 = stateOld(k, SDV_STRESS_S22) + twomu * strainInc(k,2) + alamda * trace
      s33 = stateOld(k, SDV_STRESS_S33) + twomu * strainInc(k,3) + alamda * trace
      s12 = stateOld(k, SDV_STRESS_S12) + twomu * strainInc(k,4)
      IF (nshr > 1) THEN
        s13 = stateOld(k, SDV_STRESS_S13) + twomu * strainInc(k,5)
        s23 = stateOld(k, SDV_STRESS_S23) + twomu * strainInc(k,6)
      END IF

    ! Tensão desviadora e Tensão de Von Mises
      smean = third * (s11 + s22 + s33)
      s11 = s11 - smean
      s22 = s22 - smean
      s33 = s33 - smean

      IF (nshr == 1) THEN
      vmises = SQRT(one_point_five * (s11*s11 + s22*s22 + s33*s33 + two*s12*s12))
      ELSE
      vmises = SQRT(one_point_five * (s11*s11 + s22*s22 + s33*s33 + two*s12*s12
     1   + two*s13*s13 + two*s23*s23))
      END IF

    ! Incremento de deformação plástica equivalente
      sigdif = vmises - yieldOld
      IF (sigdif > zero) THEN
      deqps = sigdif / (thremu + hard)
      END IF 
  !===============================================================
  ! --- 10. ATUALIZAÇÃO DE DANO (MECÂNICO) E TENSÃO ---
  !===============================================================
  ! --- Dano por dutilidade (Johnson-Cook) ---
      IF (vmises > 0.0) THEN
        sigmaJK = (smean * 3.0) / vmises ! Triaxialidade da tensão
      ELSE
        sigmaJK = 0.0_dp
      END IF
      epsilonJK = ABS(trace/3.0) ! Taxa de deformação (simplificado)
      epsilonf  = 0.0_dp
      IF (epsilonJK /= 0.0) THEN
        epsilonf = (D1 + D2*EXP(D3*sigmaJK)) * (1.0 + D4*LOG(epsilonJK))
        IF (ISNAN(epsilonf)) epsilonf = 0.0_dp
      END IF   
        
  ! Acumula dano mecânico
      IF ((epsilonf /= 0.0) .AND. (damage_new(nElement(k)) < dmax)) THEN
        tribological_damage = deqps / epsilonf
        damage_new(nElement(k)) = damage_new(nElement(k)) + tribological_damage

        IF (damage_new(nElement(k)) >= dmax) THEN
            stateNew(k, SDV_DELETE_FLAG) = delete_element_flag
            damage_new(nElement(k)) = dmax
            element_corrosion_status(nElement(k)) = flag_corrosion_inactive
            prop_update_trigger = flag_general_property_transfer
            CALL process_failed_element_neighbors(nElement(k), beta, OK)
            failure_occurred_this_increment = .TRUE.
        END IF
      END IF      

  ! --- Atualização da Tensão (Return Mapping com Dano) ---
      yieldNew = yieldOld + hard * deqps
      IF ((yieldNew + thremu * deqps) > 0.0_dp) THEN
          factor = yieldNew / (yieldNew + thremu * deqps)
      ELSE
          factor = 1.0_dp
      END IF

      stressNew(k,1) = (s11 * factor + smean) * (one - damage_new(nElement(k)))
      stressNew(k,2) = (s22 * factor + smean) * (one - damage_new(nElement(k)))
      stressNew(k,3) = (s33 * factor + smean) * (one - damage_new(nElement(k)))
      stressNew(k,4) = (s12 * factor) * (one - damage_new(nElement(k)))
      IF (nshr > 1) THEN
          stressNew(k,5) = (s13 * factor) * (one - damage_new(nElement(k)))
          stressNew(k,6) = (s23 * factor) * (one - damage_new(nElement(k)))
      END IF

  ! --- Remoção por excesso de deformação plástica ---
      eqps = stateOld(k, SDV_PEEQ) + deqps
      IF ((eqps > 0.137_dp) .AND. (strain_deletion_status(nElement(k)) /= flag_deletion_by_strain_disabled)) THEN
        stateNew(k, SDV_DELETE_FLAG) = delete_element_flag
        CALL process_failed_element_neighbors(nElement(k), beta, OK)
        CALL compute_nonlocal_property(ok)
        prop_update_trigger = flag_general_property_transfer
      END IF

    !===============================================================
    ! --- 11. ATUALIZAÇÃO DAS VARIÁVEIS DE ESTADO (stateNew) ---
    !===============================================================
      von_mises_stress(nElement(k)) = vmises
      damage_current(nElement(k))   = damage_new(nElement(k))

      stateNew(k, SDV_PEEQ)                 = eqps
      stateNew(k, SDV_TOTAL_DAMAGE)         = damage_new(nElement(k))
      stateNew(k, SDV_PITTING_FACTOR)       = elem_properties(nElement(k), PROP_IDX_PITTING)
      stateNew(k, SDV_CORROSION_STATUS)     = element_corrosion_status(nElement(k))
      stateNew(k, SDV_NONLOCAL_PROPERTY)    = nonlocal_property(nElement(k))
      stateNew(k, SDV_STRESS_S11)           = s11 * factor + smean
      stateNew(k, SDV_STRESS_S22)           = s22 * factor + smean
      stateNew(k, SDV_STRESS_S33)           = s33 * factor + smean
      stateNew(k, SDV_STRESS_S12)           = s12 * factor
      stateNew(k, SDV_STRESS_S13)           = s13 * factor
      stateNew(k, SDV_STRESS_S23)           = s23 * factor
      stateNew(k, SDV_VON_MISES)            = von_mises_stress(nElement(k))
      stateNew(k, SDV_DCOld)                = dCOld_k
      stateNew(k, SDV_DSCLd)                = dSCld_k
      stateNew(k, SDV_NONLOCAL_STRESS)      = non_local_stress(nElement(k))
      stateNew(k, SDV_SURFACE_FLAG)         = elem_properties(nElement(k), PROP_IDX_SURFACE_FLAG)
      stateNew(k, SDV_TOTAL_CORR_DAMAGE)    = corrosion_damage
      stateNew(k, SDV_CURRENT_TOTAL_DAMAGE) = damage_current(nElement(k))
      stateNew(k, SDV_JC_DAMAGE_INC)        = tribological_damage
      stateNew(k, SDV_TRIBO_DAMAGE_ACCUM)   = accumulated_pure_tribo
      stateNew(k, SDV_CORRODED_VOLUME)      = corroded_volume(nElement(k))      

  !===============================================================
  ! --- 12. ATUALIZAÇÃO DAS ENERGIAS ---
  !===============================================================
      IF (nshr == 1) THEN
        stressPower = half * (
     1 (stateOld(k, SDV_STRESS_S11) + stressNew(k,1)) * strainInc(k,1) + +
     2 (stateOld(k, SDV_STRESS_S22) + stressNew(k,2)) * strainInc(k,2) +
     3 (stateOld(k, SDV_STRESS_S33) + stressNew(k,3)) * strainInc(k,3) +
     4 (stateOld(k, SDV_STRESS_S12) + stressNew(k,4)) * strainInc(k,4) )
      ELSE
        stressPower = half * (
     1  (stateOld(k, SDV_STRESS_S11) + stressNew(k,1)) * strainInc(k,1) +
     2  (stateOld(k, SDV_STRESS_S22) + stressNew(k,2)) * strainInc(k,2) +
     3  (stateOld(k, SDV_STRESS_S33) + stressNew(k,3)) * strainInc(k,3) +
     4  (stateOld(k, SDV_STRESS_S12) + stressNew(k,4)) * strainInc(k,4) +
     5  (stateOld(k, SDV_STRESS_S13) + stressNew(k,5)) * strainInc(k,5) +
     6  (stateOld(k, SDV_STRESS_S23) + stressNew(k,6)) * strainInc(k,6) )
      END IF
      enerInternNew(k) = enerInternOld(k) + stressPower / density(k)

      plasticWorkInc = half * (yieldOld + yieldNew) * deqps
      enerInelasNew(k) = enerInelasOld(k) + plasticWorkInc / density(k)

        END DO ! Fim do laço principal sobre os pontos de integração
      END IF ! Fim da lógica principal do material
                 
      ! --- Liberação de Memória Alocada ---
      DEALLOCATE (damage_new)
      DEALLOCATE (relative_slip_increment)
      DEALLOCATE (contact_pressure)
      RETURN           
            
      END SUBROUTINE vumatXtrArg
!-----------------------------------------------------------------------      
!> @brief Calcula a tensão de escoamento e o encruamento isotrópico.
!!
!! Realiza uma interpolação linear na curva de tensão vs. deformação plástica equivalente,
!! fornecida na tabela de propriedades do material.
!-----------------------------------------------------------------------
      SUBROUTINE vuhard(syield, hard, eqplas, table, nvalue)
      USE SharedVariables
      USE AbaqusParameters

      ! --- Argumentos ---
      INTEGER,       INTENT(IN)  :: nvalue
      REAL(KIND=dp), INTENT(IN)  :: eqplas, table(2,*)
      REAL(KIND=dp), INTENT(OUT) :: syield, hard

      ! --- Variáveis Locais ---
      INTEGER :: k1
      REAL(kind=dp) :: eqpl1, eqpl0, deqpl, syiel0, syiel1, dsyiel
      REAL(kind=dp), PARAMETER :: zero = 0.0_dp

    ! --- Lógica Principal ---
    ! Inicializa os valores de saída assumindo extrapolação para a direita
    ! (usa o último valor da tabela se a deformação for maior que todas).
      syield = table(1, nvalue)
      hard   = zero

    ! Se a tabela tiver mais de um ponto, procura o intervalo correto.
      IF (nvalue > 1) THEN
        DO k1 = 1, nvalue - 1
          eqpl1 = table(2, k1 + 1)

          ! Se a deformação plástica está dentro do intervalo [eqpl0, eqpl1]
          IF (eqplas < eqpl1) THEN
              eqpl0  = table(2, k1)
              syiel0 = table(1, k1)
              syiel1 = table(1, k1 + 1)

              deqpl  = eqpl1 - eqpl0
              dsyiel = syiel1 - syiel0

            ! Calcula o módulo de encruamento (derivada)
            IF (deqpl > 0.0_dp) THEN
                hard = dsyiel / deqpl
            ELSE
                hard = 0.0_dp ! Evita divisão por zero
            END IF

            ! Interpolação linear para a tensão de escoamento
            syield = syiel0 + (eqplas - eqpl0) * hard
            EXIT ! Sai do laço pois o intervalo correto foi encontrado
          END IF
        END DO
      END IF
      RETURN
      END SUBROUTINE vuhard

!----------------------------------------------------------------------------------------------------------
!> @brief Ponto de entrada para controle e troca de dados com o Abaqus.
!!
!! Esta sub-rotina é chamada pelo Abaqus em diferentes pontos da análise (início, fim de incremento, etc.)
!! para permitir a inicialização de variáveis, leitura de dados e execução de rotinas de pós-processamento.
!----------------------------------------------------------------------------------------------------------   
      SUBROUTINE vexternaldb(lOp, i_Array, niArray, r_Array, nrArray)
      USE AbaqusParameters
      USE SharedVariables      

      ! --- Argumentos ---
      INTEGER,       INTENT(IN) :: lOp, niArray, nrArray, i_Array(niArray)
      REAL(KIND=dp), INTENT(IN) :: r_Array(nrArray)

      ! --- Variáveis Locais ---
      INTEGER :: k, p, l, ok, read_status, line_count, node_id
      INTEGER :: kStep, kInc, iStatus, nElems
      INTEGER, PARAMETER :: unit_num_props = 100, unit_num_normals = 101
      INTEGER, PARAMETER :: unit_num_coords = 102
      REAL(KIND=dp) :: nx, ny, nz
      CHARACTER(LEN=512) :: line_buffer  

      !=======================================================================
      ! --- INICIALIZAÇÃO DE VARIÁVEIS ---
      !=======================================================================
      kStep   = i_Array(I_INT_KSTEP)
      kInc    = i_Array(I_INT_KINC)
      iStatus = i_Array(I_INT_ISTATUS)
      nElems  = i_Array(I_INT_NTOTALELEMENTS)     

      !=======================================================================
      ! --- SELEÇÃO DO EVENTO DA ANÁLISE ---
      !=======================================================================
      SELECT CASE (lOp)     

      !-----------------------------------------------------------------------
      ! --- CASO: INÍCIO DA ANÁLISE ---
      !-----------------------------------------------------------------------
      CASE (J_INT_START_ANALYSIS)
      
      !=======================================================================
      ! 1. LEITURA DAS PROPRIEDADES DOS ELEMENTOS
      !=======================================================================
      ! 1.1 Abre arquivo e conta número de elementos
      OPEN(unit_num_props, file='E:\Pedro Bampa\corr\teste_v3\SimplesJob6.txt',
     &           status='old', action='read', iostat=read_status)           

      IF (read_status /= 0) THEN
          PRINT *, '============================================================'
          PRINT *, '>>> ERRO FATAL: Nao foi possivel abrir o arquivo Simples.txt'
          PRINT *, '============================================================'
          STOP 'Erro de abertura de arquivo de propriedades.'
      END IF

      line_count = 0
      DO
          READ(unit_num_props, '(A)', iostat=read_status) line_buffer
          IF (read_status /= 0) EXIT
          IF (TRIM(line_buffer) /= '') line_count = line_count + 1
      END DO
      REWIND(unit_num_props)
      
      ! 1.2 Aloca memória e lê propriedades
      n_elems = line_count
      PRINT *, 'Diagnostico: Arquivo de propriedades lido. N de elementos: ', n_elems
      
      CALL allocate_arrays(ok)
      IF (ok /= 1) STOP "Erro na alocação de arrays"

      DO k = 1, n_elems
          READ(unit_num_props, *) (elem_properties(k,l), l=1,14)
      END DO
      CLOSE(unit_num_props)



      ! 4. [NOVO] Ler conectividade nó-elemento
      CALL read_element_connectivity(ok)
      IF (ok /= 1) STOP "Erro na leitura da conectividade"

      !=======================================================================
      ! 2. LEITURA DOS VETORES NORMAIS
      !=======================================================================
      n_nodes = i_Array(I_INT_NTOTALNODES)
      PRINT *, 'Diagnostico: N total de nos no modelo: ', n_nodes

          ! ... ---------------------------------------------------------------!
    
    ! Alocar arrays adicionais para UMESHMOTION
      ALLOCATE(node_elements_count(n_nodes), stat=read_status)
      ALLOCATE(node_to_elem_conn(n_nodes, 20), stat=read_status) ! Assumindo máximo 20 elementos por nó
      ALLOCATE(nodal_velocities(n_nodes, 3), stat=read_status)
      
      IF (read_status /= 0) THEN
          PRINT *, '>>> ERRO FATAL: Falha ao alocar memoria para arrays UMESHMOTION'
          STOP 'Erro de alocacao de memoria.'
      END IF
      
      ! Inicializar arrays
      node_elements_count = 0
      node_to_elem_conn = 0
      nodal_velocities = 0.0_dp
      
      ! Construir mapeamento nó-elemento
      CALL build_node_to_element_connectivity(ok)
      IF (ok /= 1) STOP "Erro na construção da conectividade nó-elemento"
      
      ! Verificação - imprimir alguns valores para debug
      PRINT *, '=== Verificação Inicial ==='
      PRINT *, 'Normais para nó 1:', nodal_normals(1,:)
      PRINT *, 'Normais para nó último:', nodal_normals(n_nodes,:)
      PRINT *, 'Elementos conectados ao nó 1:', node_to_elem_conn(1,1:node_elements_count(1))
      PRINT *, 'Elementos conectados ao nó último:', node_to_elem_conn(n_nodes,1:node_elements_count(n_nodes))
      PRINT *, 'Conectividade do elemento 1:', element_nodes(1,1:num_nodes_per_element(1))
      PRINT *, 'Conectividade do elemento último:', element_nodes(n_elems,1:num_nodes_per_element(n_elems))

    ! ... ---------------------------------------------------------------!
      
      ! 2.1 Alocação do array de normais
      ALLOCATE(nodal_normals(n_nodes, 3), stat=read_status)
      IF (read_status /= 0) THEN
          PRINT *, '>>> ERRO FATAL: Falha ao alocar memoria para o array de normais.'
          STOP 'Erro de alocacao de memoria.'
      END IF
      PRINT *, 'Diagnostico: Array para normais alocado com sucesso.'

      ! 2.2 Leitura do arquivo de normais
      PRINT *, '------------------------------------------------------------'
      PRINT *, '--- Diagnostico: Iniciando leitura do arquivo de normais ---'

      OPEN(unit_num_normals, file='E:\Pedro Bampa\corr\teste_v3\normals.txt',
     &           status='old', action='read', iostat=read_status)

      IF (read_status /= 0) THEN
          PRINT *, '============================================================'
          PRINT *, '>>> ERRO FATAL: Nao foi possivel abrir o arquivo normals.txt'
          PRINT *, '>>> Verifique o caminho e a existencia do arquivo.'
          PRINT *, '============================================================'
          STOP 'Erro de abertura de arquivo de normais.'
      END IF

      nodal_normals = 0.0_dp

      DO
          READ(unit_num_normals, *, IOSTAT=read_status) node_id, nx, ny, nz
          IF (read_status /= 0) EXIT

          !PRINT *, 'Lendo normal: ID=', node_id, ', nx=', nx, ', ny=', ny, ', nz=', nz

          IF (node_id > 0 .AND. node_id <= n_nodes) THEN
              nodal_normals(node_id, 1) = nx
              nodal_normals(node_id, 2) = ny
              nodal_normals(node_id, 3) = nz
          ELSE
              PRINT *, '>>> AVISO: ID de no invalido ou fora dos limites do array. ID lido:', node_id
          END IF
      END DO
      CLOSE(unit_num_normals)
      
      PRINT *, '--- Diagnostico: Leitura do arquivo de normais concluida. ---'
      PRINT *, '------------------------------------------------------------'

      ALLOCATE(nodal_coordinates(n_nodes, 3), stat=read_status)
      IF (read_status /= 0) THEN
          PRINT *, '>>> ERRO FATAL: Falha ao alocar memoria para o array de coordenadas.'
          STOP 'Erro de alocacao de memoria.'
      END IF
      PRINT *, 'Diagnostico: Array para coordenadas alocado com sucesso.'

      ! 3.2 Leitura do arquivo de coordenadas
      OPEN(unit_num_coords, file='E:\Pedro Bampa\corr\teste_v3\nodal_coords.txt',
     &           status='old', action='read', iostat=read_status)

      IF (read_status /= 0) THEN
          PRINT *, '============================================================'
          PRINT *, '>>> ERRO FATAL: Nao foi possivel abrir o arquivo nodal_coords.txt'
          PRINT *, '============================================================'
          STOP 'Erro de abertura de arquivo de coordenadas.'
      END IF

      nodal_coordinates = 0.0_dp

      DO node_id = 1, n_nodes
          READ(unit_num_coords, *, IOSTAT=read_status) nodal_coordinates(node_id, 1), 
     &                                               nodal_coordinates(node_id, 2),
     &                                               nodal_coordinates(node_id, 3)
          IF (read_status /= 0) THEN
              PRINT *, '>>> ERRO: Leitura de coordenadas falhou no no', node_id
              EXIT
          END IF
      END DO
      CLOSE(unit_num_coords)
      
      PRINT *, 'Diagnostico: Coordenadas nodais carregadas com sucesso.'

      !=======================================================================
      ! 3. INICIALIZAÇÃO DAS VARIÁVEIS DE SIMULAÇÃO
      !=======================================================================
      ! 3.1 Configuração inicial da superfície
      CALL block_surface(ok)
      
      ! 3.2 Inicialização dos arrays temporários
      DO p = 1, n_elems
          IF (elem_properties(p,PROP_IDX_SURFACE_FLAG) == 1.0_dp) THEN
              temp_pitting_write(p) = elem_properties(p,PROP_IDX_PITTING)
              temp_surface_flag_write(p) = elem_properties(p,PROP_IDX_SURFACE_FLAG)
          ELSE
              elem_properties(p,PROP_IDX_PITTING) = 0.0_dp
              temp_pitting_write(p) = elem_properties(p,PROP_IDX_PITTING)
          END IF
      END DO
      
      ! 3.3 Configuração inicial dos cálculos não-locais
      prop_update_trigger = flag_general_property_transfer
      CALL build_influence_map(ok)
      CALL compute_nonlocal_property(ok)
      CALL transfer_properties(ok)

      !-----------------------------------------------------------------------
      ! --- CASO: FIM DE INCREMENTO ---
      !-----------------------------------------------------------------------
      CASE (J_INT_END_INCREMENT)
      ! Atualiza propriedades se houve falha de elemento
      IF (failure_occurred_this_increment) THEN
          CALL transfer_properties(ok)
          CALL compute_nonlocal_property(ok)
          failure_occurred_this_increment = .FALSE.
      END IF

      CALL compute_mass_loss(ok)
      CALL compute_nonlocal_stress(ok)

      !-----------------------------------------------------------------------
      ! --- CASO: FIM DA ANÁLISE ---
      !-----------------------------------------------------------------------
      CASE (J_INT_END_ANALYSIS)
      CALL deallocate_arrays(ok)
      CLOSE(unit_num_props)

      END SELECT
      
      RETURN
      END SUBROUTINE vexternaldb
      
!----------------------------------------------------------------------------------------------------------------------------
!> @brief Atualiza uma propriedade para os vizinhos de um elemento que falhou.
!!
!! Propaga o efeito da falha de um elemento ('nElement') para seus Vizinhos diretos. Opera em uma coluna específica da matriz
!! 'elem_properties' e armazena o resultado em um array temporário, evitando duplicação de código.
!----------------------------------------------------------------------------------------------------------------------------
      SUBROUTINE update_neighbor_property(nElement, beta, prop_column_index, target_write_array, ok)
        USE SharedVariables
        IMPLICIT NONE

      ! --- Argumentos ---
        INTEGER,       INTENT(IN)    :: nElement
        REAL(KIND=dp), INTENT(IN)    :: beta
        INTEGER,       INTENT(IN)    :: prop_column_index
        REAL(KIND=dp), INTENT(INOUT) :: target_write_array(n_elems)
        INTEGER,       INTENT(OUT)   :: ok

      ! --- Variáveis Locais ---
        INTEGER       :: j, neighbor_element
        REAL(KIND=dp) :: neighbor_prop_old, prop_from_failed

      ! --- Lógica Principal ---
      ! Pega a propriedade do elemento que falhou para propagar aos vizinhos.
        prop_from_failed = elem_properties(nElement, prop_column_index)

      ! Itera sobre os 6 vizinhos do elemento que falhou (pressupõe malha estruturada).
        DO j = 1, max_face_neighbors
          ! As colunas de 2 a 7 de elem_properties contêm os IDs dos vizinhos.
          neighbor_element = elem_properties(nElement, j + 1)

          ! Procede apenas se o vizinho existir (ID diferente de zero).
          IF (neighbor_element /= 0) THEN
              ! Se o vizinho não estiver "travado", calcula e atualiza a propriedade.
            IF (property_update_lock(neighbor_element) /= key_prop_update_locked) THEN
                neighbor_prop_old = elem_properties(neighbor_element, prop_column_index)

                ! A lógica MAX garante que a propriedade do vizinho apenas aumente ou
                ! permaneça a mesma, nunca diminua com a falha de outro vizinho.
                target_write_array(neighbor_element) = MAX(neighbor_prop_old, prop_from_failed * beta)
            END IF
          END IF
        END DO

      ! Zera a propriedade no elemento que efetivamente falhou e o trava para
      ! futuras atualizações de vizinhança.
        target_write_array(nElement) = 0.0_dp
        property_update_lock(nElement) = key_prop_update_locked

        ok = 1
        RETURN
      END SUBROUTINE update_neighbor_property
!-----------------------------------------------------------------------------------------------------
!> @brief Orquestra todas as atualizações para os vizinhos de um elemento que falhou.
!!
!! Quando um elemento falha na VUMAT, esta rotina é chamada para garantir que as consequências da falha
!! (atualização de pitting e flag de superfície) sejam propagadas corretamente aos vizinhos.
!-----------------------------------------------------------------------------------------------------
      SUBROUTINE process_failed_element_neighbors(nElement, beta, ok)
        USE SharedVariables
        IMPLICIT NONE

        ! --- Argumentos ---
        INTEGER,       INTENT(IN)  :: nElement
        REAL(KIND=dp), INTENT(IN)  :: beta
        INTEGER,       INTENT(OUT) :: ok
        
        ! --- Lógica ---
        ! Esta rotina encapsula as duas ações que sempre ocorrem juntas na falha.

        ! 1. Atualiza a propriedade de flag de superfície (coluna 8).
        CALL update_neighbor_property(nElement, beta, PROP_IDX_SURFACE_FLAG, temp_surface_flag_write, ok)

        ! 2. Atualiza a propriedade de fator de pitting (coluna 9).
        CALL update_neighbor_property(nElement, beta, PROP_IDX_PITTING, temp_pitting_write, ok)
        
        ! Retorna OK. A lógica de tratamento de erro pode ser aprimorada.
        ok = 1  
        RETURN
      END SUBROUTINE process_failed_element_neighbors
! -------------------------------------------------------------------------
!> @brief Calcula a propriedade não-local (homogeneizada) para elementos de superfície.
!!
!! Implementa um modelo de média ponderada para calcular uma propriedade não-local,
!! considerando a influência dos elementos vizinhos dentro de um raio 'lr'.
! -------------------------------------------------------------------------
      SUBROUTINE compute_nonlocal_property(ok)
        USE SharedVariables
        IMPLICIT NONE

        ! --- Argumentos ---
        INTEGER, INTENT(OUT) :: ok

        ! --- Variáveis Locais ---
        INTEGER       :: i, p, neighbor_element
        REAL(KIND=dp) :: denominator_sum, numerator_sum, weight_p
        REAL(KIND=dp) :: dist_sq, lr_sq

        ! --- Lógica Principal ---
        ! Pré-calcula o quadrado do raio de influência para otimização.
        lr_sq = lr**2.0_dp

        DO i = 1, n_elems
            ! Procede somente se o elemento 'i' for de superfície (propriedade na coluna 9 não é nula).
            IF (elem_properties(i, PROP_IDX_PITTING) /= 0.0_dp) THEN

            ! Inicializa as somas com a contribuição do próprio elemento 'i'.
              denominator_sum = elem_properties(i, PROP_IDX_VOLUME) ! Denominador: Soma ponderada dos volumes
              numerator_sum   = elem_properties(i, PROP_IDX_PITTING) * elem_properties(i, PROP_IDX_VOLUME) ! Numerador: Soma ponderada de (prop * vol)

            ! Itera sobre os vizinhos para calcular suas contribuições.
              DO p = 1, counter(i)
              neighbor_element = influ(i, p)

              ! Considera apenas vizinhos que também são de superfície.
                IF (elem_properties(neighbor_element, PROP_IDX_PITTING) /= 0.0_dp) THEN
                    dist_sq = dist(i, p)**2.0_dp
                        
                    ! Função de peso quadrática
                    weight_p = (1.0_dp - (dist_sq / lr_sq))**2.0_dp

                    denominator_sum = denominator_sum + (weight_p * elem_properties(neighbor_element, PROP_IDX_VOLUME))
                    numerator_sum = numerator_sum + (weight_p * elem_properties(neighbor_element, PROP_IDX_PITTING)
     1                                                        * elem_properties(neighbor_element, PROP_IDX_VOLUME))
                END IF
              END DO

                ! Calcula o valor final da propriedade não-local (média ponderada).
                IF (denominator_sum > 1.0e-20_dp) THEN
                    temp_nonlocal_write(i) = numerator_sum / denominator_sum
                ELSE
                    temp_nonlocal_write(i) = 0.0_dp
                END IF
            ELSE
                ! Para elementos internos, a propriedade não-local é zero.
                temp_nonlocal_write(i) = 0.0_dp
            END IF
        END DO

        ok = 1
        RETURN
      END SUBROUTINE compute_nonlocal_property
!-----------------------------------------------------------------------------------------------------------------
!> !> @brief Transfere propriedades calculadas dos arrays temporários para os globais.
!!
!! Esta rotina atua como um portão, copiando dados dos vetores temporários (ex: temp_nonlocal_write)
!! para as variáveis de estado permanentes (ex: nonlocal_property). A execução é controlada por um gatilho global.
!-----------------------------------------------------------------------------------------------------------------
      SUBROUTINE transfer_properties(ok)
          USE SharedVariables
          IMPLICIT NONE

          ! --- Argumentos ---
          INTEGER, INTENT(OUT) :: ok

          ! --- Variáveis Locais ---
          INTEGER :: i

          ! --- Lógica Principal ---
          ! A transferência só ocorre se o gatilho global estiver ativado.
          IF (prop_update_trigger == flag_general_property_transfer) THEN

              ! Copia os dados dos arrays temporários para os arrays globais.
              DO i = 1, n_elems
                  nonlocal_property(i) = temp_nonlocal_write(i)
                  elem_properties(i,PROP_IDX_PITTING) = temp_pitting_write(i)
                  elem_properties(i,PROP_IDX_SURFACE_FLAG) = temp_surface_flag_write(i)
              END DO

              ! Desativa o gatilho para evitar reexecução desnecessária.
              ! Ele será reativado quando uma nova falha de elemento ocorrer.
              prop_update_trigger = flag_general_property_transfer_locked
          END IF

          ok = 1
          RETURN
      END SUBROUTINE transfer_properties   
!------------------------------------------------------------------------------------------------------------------
!> @brief Calcula o mapa de influência para cada elemento do modelo.
!!
!! Para cada elemento 'i', esta rotina encontra todos os outros elementos 'j'cuja distância euclidiana
!! seja menor que o comprimento intrínseco 'lr'. Ela armazena os IDs, as distâncias e a contagem total de vizinhos.
!------------------------------------------------------------------------------------------------------------------    
      SUBROUTINE build_influence_map(ok)
        USE SharedVariables
        IMPLICIT NONE

        ! --- Argumentos ---
        INTEGER, INTENT(OUT) :: ok

        ! --- Variáveis Locais ---
        INTEGER       :: t, u, i, j, m
        REAL(KIND=dp) :: distance

        ! --- 1. Inicializa os vetores e matrizes de influência ---
        DO t = 1, n_elems
            counter(t) = 0
            DO u = 1, max_influence_elems
                dist(t,u)  = 0.0_dp
                influ(t,u) = 0
            END DO
        END DO

        ! --- 2. Calcula a distância entre cada par de elementos (i, j) ---
        ! Este é um laço O(N^2), pode ser lento para malhas muito grandes.
        DO i = 1, n_elems
            m = 1 ! Reinicia o contador de vizinhos para o elemento 'i'
            DO j = 1, n_elems
                IF (j /= i) THEN

                ! Cálculo da distância euclidiana
                distance = SQRT((elem_properties(j,PROP_IDX_COORD_X) - elem_properties(i,PROP_IDX_COORD_X))**2.0 +
     1                          (elem_properties(j,PROP_IDX_COORD_Y) - elem_properties(i, PROP_IDX_COORD_Y))**2.0 +
     2                          (elem_properties(j,PROP_IDX_COORD_Z) - elem_properties(i,PROP_IDX_COORD_Z))**2.0 )

                ! Se 'j' estiver dentro do raio de influência, armazena-o como vizinho.
                  IF (distance < lr) THEN
                          ! Verifica se o número de vizinhos não excede o limite alocado.
                          IF (m <= max_influence_elems) THEN
                              influ(i,m)   = j
                              dist(i,m)    = distance
                              counter(i)   = m
                              m = m + 1
                          ELSE                        
                          END IF
                  END IF
                END IF
            END DO
        END DO

        ok = 1
        RETURN
      END SUBROUTINE build_influence_map

!--------------------------------------------------------------------------------------------------------------
!> @brief Reserva memória para as variáveis compartilhadas de forma segura.
!!
!! Garante que a memória só seja alocada se ainda não estiver, evitando
!! erros de "array is already allocated".
!--------------------------------------------------------------------------------------------------------------
      SUBROUTINE allocate_arrays(ok)
          USE SharedVariables
          IMPLICIT NONE

          ! --- Argumentos ---
          INTEGER, INTENT(OUT) :: ok

          ! --- Variáveis Locais ---
          INTEGER :: stat ! Variável para verificar o status da alocação

          ! --- Alocação Segura de todos os Arrays Globais ---
          IF (.NOT. ALLOCATED(elem_properties))          ALLOCATE(elem_properties(n_elems, 14), stat=stat)
          IF (.NOT. ALLOCATED(nonlocal_property))        ALLOCATE(nonlocal_property(n_elems), stat=stat)
          IF (.NOT. ALLOCATED(property_update_lock))     ALLOCATE(property_update_lock(n_elems), stat=stat)
          IF (.NOT. ALLOCATED(influ))                    ALLOCATE(influ(n_elems, max_influence_elems), stat=stat)
          IF (.NOT. ALLOCATED(dist))                     ALLOCATE(dist(n_elems, max_influence_elems), stat=stat)
          IF (.NOT. ALLOCATED(temp_pitting_write))       ALLOCATE(temp_pitting_write(n_elems), stat=stat)
          IF (.NOT. ALLOCATED(temp_surface_flag_write))  ALLOCATE(temp_surface_flag_write(n_elems), stat=stat)
          IF (.NOT. ALLOCATED(counter))                  ALLOCATE(counter(n_elems), stat=stat)
          IF (.NOT. ALLOCATED(temp_nonlocal_write))      ALLOCATE(temp_nonlocal_write(n_elems), stat=stat)
          IF (.NOT. ALLOCATED(damage_current))           ALLOCATE(damage_current(n_elems), stat=stat)
          IF (.NOT. ALLOCATED(element_corrosion_status)) ALLOCATE(element_corrosion_status(n_elems), stat=stat)
          IF (.NOT. ALLOCATED(strain_deletion_status))   ALLOCATE(strain_deletion_status(n_elems), stat=stat)
          IF (.NOT. ALLOCATED(von_mises_stress))         ALLOCATE(von_mises_stress(n_elems), stat=stat)
          IF (.NOT. ALLOCATED(non_local_stress))         ALLOCATE(non_local_stress(n_elems), stat=stat)
          IF (.NOT. ALLOCATED(element_recession_velocity)) ALLOCATE(element_recession_velocity(n_elems), stat=stat)
          IF (.NOT. ALLOCATED(initial_nodal_normals))      ALLOCATE(initial_nodal_normals(n_elems, 3), stat=stat)
          IF (.NOT. ALLOCATED(element_nodes))              ALLOCATE(element_nodes(n_elems, max_nodes_per_elem), stat=stat)
          IF (.NOT. ALLOCATED(num_nodes_per_element))      ALLOCATE(num_nodes_per_element(n_elems), stat=stat)

          ! Verificar o status da alocação
           IF (stat /= 0) THEN
               PRINT *, "Erro ao alocar memoria!"
               STOP
           END IF

          ok = 1
          RETURN
      END SUBROUTINE allocate_arrays

!--------------------------------------------------------------------------------------------------------------
!> @brief Desaloca memória das variáveis compartilhadas ao fim da análise.
!!
!! Garante a liberação de toda a memória alocada dinamicamente para evitar vazamentos de memória (memory leaks).
!--------------------------------------------------------------------------------------------------------------           
      SUBROUTINE deallocate_arrays(ok)
        USE SharedVariables
        IMPLICIT NONE
        INTEGER, intent(out) :: ok

        ! --- Desalocação Segura de todos os Arrays Globais ---
        IF (allocated(elem_properties))          DEALLOCATE(elem_properties)
        IF (allocated(nonlocal_property))        DEALLOCATE(nonlocal_property)
        IF (allocated(property_update_lock))     DEALLOCATE(property_update_lock)
        IF (allocated(influ))                    DEALLOCATE(influ) 
        IF (allocated(dist))                     DEALLOCATE(dist)
        IF (allocated(temp_pitting_write))       DEALLOCATE(temp_pitting_write)
        IF (allocated(temp_surface_flag_write))  DEALLOCATE(temp_surface_flag_write)
        IF (allocated(counter))                  DEALLOCATE(counter)
        IF (allocated(temp_nonlocal_write))      DEALLOCATE(temp_nonlocal_write)
        IF (allocated(damage_current))           DEALLOCATE(damage_current)
        IF (allocated(element_corrosion_status)) DEALLOCATE(element_corrosion_status)
        IF (allocated(strain_deletion_status))   DEALLOCATE(strain_deletion_status)
        IF (allocated(von_mises_stress))         DEALLOCATE(von_mises_stress)
        IF (allocated(non_local_stress))         DEALLOCATE(non_local_stress)
        IF (allocated(element_recession_velocity)) DEALLOCATE(element_recession_velocity)
        IF (allocated(initial_nodal_normals))      DEALLOCATE(initial_nodal_normals)
        IF (allocated(element_nodes))              DEALLOCATE(element_nodes)
        IF (allocated(num_nodes_per_element))      DEALLOCATE(num_nodes_per_element)

        ok=1
      END SUBROUTINE deallocate_arrays    
!---------------------------------------------------------------------------------------------------------------------------------
!> @brief Calcula a massa total perdida por corrosão a cada incremento.
!!
!! Soma a perda de massa em todos os elementos com base no dano atual e imprime o total no console se a mudança for significativa.
!---------------------------------------------------------------------------------------------------------------------------------
      SUBROUTINE compute_mass_loss(ok)
        USE SharedVariables
        IMPLICIT NONE

        ! --- Argumentos ---
        INTEGER, INTENT(OUT) :: ok

        ! --- Variáveis Locais ---
        REAL(KIND=dp) :: rho, deltam, dmCO
        INTEGER       :: i

        ! --- Parâmetros Locais ---
        ! Nota: A densidade (rho) é definida localmente aqui.
        rho    = 7.86e-09 ! [ton/mm^3] Densidade do material
        deltam = 0.0_dp   ! Limiar de mudança de massa para impressão

        ! --- Lógica Principal ---
        ! Zera a massa corroída do incremento atual antes de somar.
        mCO = 0.0_dp

        ! Soma a contribuição de massa perdida de cada elemento.
        ! (Ver sugestões sobre o índice 10)
        DO i = 1, n_elems
            mCO = mCO + (damage_current(i) * rho * elem_properties(i,PROP_IDX_VOLUME))
        END DO

        ! Calcula a variação de massa no incremento.
        dmCO = mCO - lmCO

        ! Se a corrosão não foi parada e a mudança de massa é significativa,
        ! imprime o status e atualiza a massa de referência.
        IF((global_corrosion_stop_flag /= flag_global_stop) .AND. (dmCO >= deltam)) THEN
            !WRITE(*,*) 'Massa Corroída Total:', mCO, ' Tempo:', timetot
            lmCO = mCO ! Atualiza o valor de massa do passo anterior.
        END IF

        ok = 1
        RETURN
      END SUBROUTINE compute_mass_loss
!--------------------------------------------------------------------------------------------------------
!> @brief Bloqueia a corrosão em regiões específicas da geometria.
!!
!! Chamada no início da análise para inicializar os flags de corrosão com base na posição Z dos elementos.
!--------------------------------------------------------------------------------------------------------
      SUBROUTINE block_surface(ok)
          USE SharedVariables
          IMPLICIT NONE

          ! --- Argumentos ---
          INTEGER, INTENT(OUT) :: ok

          ! --- Variáveis Locais ---
          INTEGER       :: k
          REAL(KIND=dp) :: coordz1, coordz2

          ! --- Lógica Principal ---
          ! Controle geométrico para definir a zona de não-corrosão.
          ! (Ver sugestões sobre os valores e índices)
          coordz1 = 16.8
          coordz2 = -0.2

          DO k = 1, n_elems
              ! Se o centroide do elemento está fora da região de interesse...
              IF ((elem_properties(k,13) >= coordz1) .OR. (elem_properties(k,13) <= coordz2)) THEN
                  ! ...desativa todas as flags relevantes para corrosão e dano.
                  property_update_lock(k)                         = key_prop_update_locked
                  elem_properties(k,PROP_IDX_SURFACE_FLAG)        = 0.0_dp
                  element_corrosion_status(k)                     = flag_corrosion_inactive
                  damage_current(k)                               = 0.0_dp
                  strain_deletion_status(k)                       = flag_deletion_by_strain_disabled
              ELSE
                  ! ...senão, ativa as flags para permitir corrosão e dano.
                  property_update_lock(k)     = key_prop_update_unlocked
                  element_corrosion_status(k) = flag_corrosion_active
                  damage_current(k)           = 0.0_dp
                  strain_deletion_status(k)   = flag_deletion_by_strain_enabled
              END IF
          END DO

          ok = 1
          RETURN
      END SUBROUTINE block_surface
!--------------------------------------------------------------------------------------------------------------------------
!> @brief Interrompe globalmente o processo de corrosão.
!!
!! Chamada quando uma condição de parada é atingida (ex: massa corroída máxima). Desativa a corrosão para TODOS os elementos.
!--------------------------------------------------------------------------------------------------------------------------
      SUBROUTINE block_all(ok)
        USE SharedVariables
        IMPLICIT NONE

        ! --- Argumentos ---
        INTEGER, INTENT(OUT) :: ok

        ! --- Variáveis Locais ---
        INTEGER :: k

        ! --- Lógica Principal ---
        DO k = 1, n_elems
            ! Trava a atualização de propriedades de vizinhança.
            property_update_lock(k) = key_prop_update_locked
            ! Desativa o cálculo de corrosão para o elemento.
            element_corrosion_status(k) = flag_corrosion_inactive
        END DO

        ! Ativa o flag global para sinalizar a parada completa da corrosão.
        global_corrosion_stop_flag = flag_global_stop

        ok = 1
        RETURN
      END SUBROUTINE block_all
!-----------------------------------------------------------------------
!> @brief Calcula a tensão não-local para os elementos de superfície.
!!
!! Implementa uma média ponderada não-local para a tensão, análoga à
!! sub-rotina `compute_nonlocal_property`. O resultado só é armazenado se
!! exceder a tensão de escoamento do material (`sigth`).
!-----------------------------------------------------------------------
      SUBROUTINE compute_nonlocal_stress(ok)
          USE SharedVariables
          IMPLICIT NONE

          ! --- Argumentos ---
          INTEGER, INTENT(OUT) :: ok

          ! --- Variáveis Locais ---
          INTEGER       :: i, p, neighbor_element
          REAL(KIND=dp) :: denominator_sum, numerator_sum, weight_p
          REAL(KIND=dp) :: dist_sq, lr_sq

          ! --- Lógica Principal ---
          lr_sq = lr**2.0_dp

          DO i = 1, n_elems
              ! Considera apenas elementos de superfície não bloqueados.
            IF ((elem_properties(i, PROP_IDX_PITTING) /= 0.0_dp) .AND.
     1                (element_corrosion_status(i) /= flag_corrosion_inactive)) THEN             

                  ! Inicializa as somas com a contribuição do próprio elemento 'i'.
                  denominator_sum = elem_properties(i, PROP_IDX_VOLUME)
                  numerator_sum   = von_mises_stress(i) * elem_properties(i, PROP_IDX_VOLUME)

                  ! Itera sobre os vizinhos para calcular suas contribuições.
                  DO p = 1, counter(i)
                      neighbor_element = influ(i, p)

                      IF (element_corrosion_status(neighbor_element) /= flag_corrosion_inactive) THEN
                          dist_sq  = dist(i, p)**2.0_dp
                          weight_p = (1.0_dp - (dist_sq / lr_sq))**2.0_dp

                          denominator_sum = denominator_sum + (weight_p * elem_properties(neighbor_element, PROP_IDX_VOLUME))
                          numerator_sum = numerator_sum + (weight_p * von_mises_stress(neighbor_element) *
     1                elem_properties(neighbor_element, PROP_IDX_VOLUME))
                      END IF
                  END DO

                  ! Calcula o valor final da tensão não-local (média ponderada).
                  IF (denominator_sum > 1.0e-20_dp) THEN
                      non_local_stress(i) = numerator_sum / denominator_sum
                  ELSE
                      non_local_stress(i) = 0.0_dp
                  END IF

                  ! Aplica o limiar da tensão de escoamento (threshold).
                  IF (non_local_stress(i) <= sigth) THEN
                      non_local_stress(i) = 0.0_dp
                  END IF
              ELSE
                  ! Para elementos internos ou inativos, a tensão não-local é zero.
                  non_local_stress(i) = 0.0_dp
              END IF
          END DO

          ok = 1
          RETURN
      END SUBROUTINE compute_nonlocal_stress

!-----------------------------------------------------------------------
!> @brief 
!!
!! 
!! 
!! 
!-----------------------------------------------------------------------


      SUBROUTINE read_element_connectivity(ok)
          USE SharedVariables
          IMPLICIT NONE
          INTEGER, INTENT(OUT) :: ok
          INTEGER :: unit_num, i, elem_id, stat
          CHARACTER(LEN=256) :: filename
          
          filename = 'E:\Pedro Bampa\corr\teste_v3\element_connectivity.txt'  ! Arquivo pré-gerado
          unit_num = 50
           ! Debug: Início da subrotina
           WRITE(6, *) "=== read_element_connectivity chamada ==="

          OPEN(unit_num, FILE=filename, STATUS='OLD', ACTION='READ', IOSTAT=stat)
          IF (stat /= 0) THEN
              PRINT *, "Erro abrindo arquivo de conectividade"
              ok = 0
              RETURN
          END IF
          
          DO i = 1, n_elems
              READ(unit_num, *, IOSTAT=stat) elem_id, element_nodes(elem_id,:)
              IF (stat /= 0) EXIT
              num_nodes_per_element(elem_id) = COUNT(element_nodes(elem_id,:) > 0)
          END DO
          
          CLOSE(unit_num)
          ok = 1
      END SUBROUTINE read_element_connectivity
!-----------------------------------------------------------------------
!> @brief 
!!
!! 
!! 
!! 
!-----------------------------------------------------------------------
      SUBROUTINE build_node_to_element_connectivity(ok)
      USE SharedVariables
      IMPLICIT NONE
      INTEGER, INTENT(OUT) :: ok
      INTEGER :: elem, node_id, i, stat
      REAL(KIND=dp) :: debug_coords(3)
      
      ! Primeiro: determinar o máximo de elementos por nó
      INTEGER :: max_elems_per_node = 0
      INTEGER, ALLOCATABLE :: temp_count(:)

      ! Variáveis para debug
      INTEGER, SAVE :: call_count = 0
      CHARACTER(LEN=100) :: debug_message
      
      ! Incrementa o contador de chamadas
      call_count = call_count + 1
      
      ! Mensagem de debug
      debug_coords = 0.0_dp  ! Ou valores reais se disponíveis
      WRITE(debug_message, '(A, I6, A, I6, A, 3F10.5)') 
     &     'UMESHMOTION CALL #', call_count,
     &     ' | Node ID: ', FIND_CLOSEST_NODE(debug_coords),
     &     ' | Coords: ', debug_coords(1), debug_coords(2), debug_coords(3)
      
      ! Escreve no arquivo .msg (usando unit=6, que é o output padrão do Abaqus)
      WRITE(6, *) trim(debug_message)
      
      ALLOCATE(temp_count(n_nodes), stat=stat)
      temp_count = 0
      
      ! Passada 1: contar elementos por nó
      DO elem = 1, n_elems
          DO i = 1, num_nodes_per_element(elem)
              node_id = element_nodes(elem, i)
              IF (node_id > 0 .AND. node_id <= n_nodes) THEN
                  temp_count(node_id) = temp_count(node_id) + 1
              END IF
          END DO
      END DO
      
      max_elems_per_node = MAXVAL(temp_count)
      
      ! Redimensionar arrays se necessário
      IF (ALLOCATED(node_to_elem_conn)) THEN
          IF (SIZE(node_to_elem_conn, 2) < max_elems_per_node) THEN
              DEALLOCATE(node_to_elem_conn)
              ALLOCATE(node_to_elem_conn(n_nodes, max_elems_per_node), stat=stat)
          END IF
      ELSE
          ALLOCATE(node_to_elem_conn(n_nodes, max_elems_per_node), stat=stat)
      END IF
      
      ! Passada 2: preencher conectividade
      node_elements_count = 0
      DO elem = 1, n_elems
          DO i = 1, num_nodes_per_element(elem)
              node_id = element_nodes(elem, i)
              IF (node_id > 0 .AND. node_id <= n_nodes) THEN
                  node_elements_count(node_id) = node_elements_count(node_id) + 1
                  IF (node_elements_count(node_id) <= max_elems_per_node) THEN
                      node_to_elem_conn(node_id, node_elements_count(node_id)) = elem
                  END IF
              END IF
          END DO
      END DO
      
      DEALLOCATE(temp_count)
      ok = 1
      
      ! Verificação
      PRINT *, 'Conectividade construída. Máx elementos/nó:', max_elems_per_node
      END SUBROUTINE build_node_to_element_connectivity
!-----------------------------------------------------------------------
!> @brief 
!!
!! 
!! 
!! 
!-----------------------------------------------------------------------
              SUBROUTINE UMESHMOTION(
     1     U, DU, V, A, JTYPE, TIME, DTIME, KSTEP, KINC,
     2     JELEM, LAYERP, KSPT, COORDS, JLENGTH, NNODE, JDOF,
     3     NDI, NSHR, KFIELD, KPROPS, KJPROPS, KMAT, KDLOAD,
     4     KPREDEF, NPREDF, JPROPS, JDLTYP, ADLMAG, PREDEF,
     5     DPROPS, TDELT, DDLMAG, MDLOAD, PNEWDT, JPRED,
     6     LFLAGS, DVRESULT)
     
      USE SharedVariables
      USE AbaqusParameters
      !IMPLICIT NONE
      
      ! Declaração de argumentos
      DIMENSION U(3), DU(3,3), V(3), A(3)
      DIMENSION JDLTYP(MDLOAD,*), ADLMAG(MDLOAD,*)
      DIMENSION DDLMAG(MDLOAD,*), PREDEF(2,NPREDF,*)
      DIMENSION DPROPS(*), JPROPS(*), DVRESULT(7)
      DIMENSION LFLAGS(*)
      REAL(KIND=dp) :: coords_dp(3)
      
      ! Variáveis locais
      INTEGER :: i, elem, node_id
      REAL(KIND=dp) :: node_velocity(3), weight, total_weight
      REAL(KIND=dp) :: min_dist, distance, node_coords(3)
      LOGICAL :: node_found
      
      ! Inicialização
      DVRESULT = 0.0_dp
      node_id = 0  ! Valor padrão para caso de erro
      node_found = .FALSE.
      
      IF (ALLOCATED(nodal_coordinates)) THEN
          ! Chama a função que agora está dentro do módulo SharedVariables
          coords_dp = REAL(COORDS, KIND=dp) 
          node_id = FIND_CLOSEST_NODE(coords_dp)
          
          ! Verificação adicional de segurança
          IF (node_id > 0 .AND. node_id <= SIZE(nodal_coordinates,1)) THEN
              node_found = .TRUE.
          ELSE
              WRITE(*,*) 'AVISO: Nó não encontrado para coordenadas:', COORDS
          END IF
      ELSE
          WRITE(*,*) 'ERRO: nodal_coordinates não alocado em UMESHMOTION'
      END IF

! 2. Continuar apenas
      IF (node_found) THEN 
          total_weight = 0.0_dp
          node_velocity = 0.0_dp
          
          DO i = 1, node_elements_count(node_id)
              elem = node_to_elem_conn(node_id, i)
              weight = SQRT(SUM(nodal_normals(node_id,:)**2))
              total_weight = total_weight + weight
              
              node_velocity = node_velocity + element_recession_velocity(elem) *
     1              nodal_normals(node_id,:) * weight
          END DO
          
          IF (total_weight > 0.0_dp) THEN
              node_velocity = node_velocity / total_weight
          END IF
          
          nodal_velocities(node_id,:) = node_velocity
          DVRESULT(1:3) = node_velocity
          WRITE(*,*) 'UMESHMOTION called for node:', node_id, 'Velocity:', DVRESULT(1:3)
      END IF
      
      RETURN
      END SUBROUTINE UMESHMOTION



        