
    module Reionization
    use Precision
    use MiscUtils
    use classes
    use results
    implicit none
    private
    
    !Default tanh reionization, and an alternative exponential model with fixed minimum z_re
    
    !This module has smooth tanh reionization of specified mid-point (z_{re}) and width
    !The tanh function is in the variable (1+z)**Rionization_zexp
    !Rionization_zexp=1.5 has the property that for the same z_{re}
    !the optical depth agrees with infinitely sharp model for matter domination
    !So tau and zre can be mapped into each other easily (for any symmetric window)
    !However for generality the module maps tau into z_{re} using a binary search
    !so could be easily modified for other monatonic parameterizations.
    
    !The ionization history must be twice differentiable.

    !AL March 2008
    !AL July 2008 - added trap for setting optical depth without use_optical_depth
    !AL Aug 2023 - added exponential model and refactored classes
    
    !See CAMB notes for further discussion: http://cosmologist.info/notes/CAMB.pdf

    real(dl), parameter :: Reionization_DefFraction = -1._dl
    !if -1 set from YHe assuming Hydrogen and first ionization of Helium follow each other

    real(dl) :: Tanh_zexp = 1.5_dl

    type, extends(TReionizationModel) :: TBaseTauWithHeReionization
        ! Parameterization that can take tau as an input, using redshift as a one-parameter mapping to tau
        ! includes simple tanh fitting of second reionization of helium
        logical    :: use_optical_depth = .false.
        real(dl)   :: redshift = 10._dl
        real(dl)   :: optical_depth = 0._dl
        real(dl)   :: fraction = Reionization_DefFraction
        !Parameters for the second reionization of Helium
        logical    :: include_helium_fullreion  = .true.
        real(dl)   :: helium_redshift  = 3.5_dl
        real(dl)   :: helium_delta_redshift  = 0.4_dl
        real(dl)   :: helium_redshiftstart  = 5.5_dl
        real(dl)   :: tau_solve_accuracy_boost = 1._dl
        real(dl)   :: timestep_boost =  1._dl
        real(dl)   :: max_redshift = 50._dl
        real(dl)   :: min_redshift = 0._dl
        !The rest are internal to this module.
        real(dl), private ::  fHe
        class(CAMBdata), pointer :: State
    contains
    procedure :: ReadParams => TBaseTauWithHeReionization_ReadParams
    procedure :: Init => TBaseTauWithHeReionization_Init
    procedure, nopass ::  GetZreFromTau => TBaseTauWithHeReionization_GetZreFromTau
    procedure, private :: zreFromOptDepth => TBaseTauWithHeReionization_zreFromOptDepth
    procedure :: SecondHelium_xe => TBaseTauWithHeReionization_SecondHelium_xe
    procedure :: SetParamsForZre => TBaseTauWithHeReionization_SetParamsForZre
    procedure :: Validate => TBaseTauWithHeReionization_Validate
    end type TBaseTauWithHeReionization

    type, extends(TBaseTauWithHeReionization) :: TTanhReionization
        real(dl)   :: delta_redshift = 0.5_dl
        !The rest are internal to this module.
        real(dl), private ::  WindowVarMid, WindowVarDelta
    contains
    procedure :: x_e => TTanhReionization_xe
    procedure :: get_timesteps => TTanhReionization_get_timesteps
    procedure :: ReadParams => TTanhReionization_ReadParams
    procedure :: Validate => TTanhReionization_Validate
    procedure :: SetParamsForZre => TTanhReionization_SetParamsForZre
    procedure, nopass :: SelfPointer => TTanhReionization_SelfPointer
    end type TTanhReionization

    type, extends(TBaseTauWithHeReionization) :: TExpReionization
        ! An ionization fraction that decreases exponentially at high z, saturating to fully inionized at fixed redshift.
        ! This model has a minimum non-zero tau
        ! Similar to e.g.  arXiv:1509.02785, arXiv:2006.16828
        real(dl)   :: reion_redshift_complete = 6.1_dl
        real(dl)   :: reion_exp_smooth_width = 0.02_dl !modifies expential at reion_redshift_complete so derivatives continuous
        real(dl)   :: reion_exp_power = 1._dl  !scaling propto exp(-lambda (z-reion_redshift_complete)**reion_exp_power) at high z
    contains
    procedure :: x_e => TExpReionization_xe
    procedure :: get_timesteps => TExpReionization_get_timesteps
    procedure :: Init => TExpReionization_Init
    procedure :: ReadParams => TExpReionization_ReadParams
    procedure, nopass :: SelfPointer => TExpReionization_SelfPointer
    end type TExpReionization

    !MT
    type, extends(TBaseTauWithHeReionization) :: TPCAReionization
        ! An ionization fraction from Principal Component Analysis
        ! This model fit amplitude of eigenvectors of xe(z)
        ! See e.g.  arXiv:0303400, arXiv:1609.04788
        real(dl) :: reion_pca_mode1=0._dl !amplitudes of PCA modes
        real(dl) :: reion_pca_mode2=0._dl !amplitudes of PCA modes
        real(dl) :: reion_pca_mode3=0._dl !amplitudes of PCA modes
        real(dl) :: reion_pca_mode4=0._dl !amplitudes of PCA modes
        real(dl) :: reion_pca_mode5=0._dl !amplitudes of PCA modes
        real(dl) :: reion_pca_mode6=0._dl !amplitudes of PCA modes
        real(dl) :: reion_pca_mode7=0._dl !amplitudes of PCA modes
        real(dl) :: reion_pca_mode8=0._dl !amplitudes of PCA modes
        real(dl) :: reion_pca_mode9=0._dl !amplitudes of PCA modes
        real(dl) :: reion_pca_mode10=0._dl !amplitudes of PCA modes
        real, dimension(40) :: reion_redshift = [5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 14.5, 15.0, 15.5, 16.0, 16.5, 17.0, 17.5, 18.0, 18.5, 19.0, 19.5, 20.0, 20.5, 21.0, 21.5, 22.0, 22.5, 23.0, 23.5, 24.0, 24.5, 25.0]
        real, dimension(400) :: reion_pca = [&
             0.01015315327687579, 0.014044269873504947, 0.0271696164906656, 0.0468109891410037, 0.07391725998892867, 0.107028188927783, 0.14667064632324023, 0.19258418260856175, 0.23820069567523205, 0.2921213328684633, 0.35417203126959773, 0.40934739981594465, 0.4681097148238654, 0.5362119235492643, 0.601855620854531, 0.6605076245382242, 0.733762741639875, 0.7901406689587728, 0.8529319596730001, 0.9150183390706853, 0.9663065937468368, 1.0197286831007244, 1.0704050448924656, 1.1115053115335163, 1.1601023572137361, 1.1900009402751504, 1.2435637163608855, 1.2691507539350393, 1.2989388240710134, 1.354701441893458, 1.3727499633831253, 1.387179477146869, 1.411186619771781, 1.4405480740231238, 1.4504854928613082, 1.489381248528841, 1.5054861946651645, 1.5422126725742036, 1.5466451659610203, 1.5564429621855236, &
             -0.8827328536970834, -0.9714199210572283, -1.1077739620331795, -1.2094811885145722, -1.319351645601336, -1.3917184085252867, -1.4716964163700572, -1.5261983649517896, -1.537673295732187, -1.5520540008634602, -1.554389204140202, -1.5150294828312023, -1.4521012393038966, -1.420068416834558, -1.3399303197776187, -1.2386156051961263, -1.1637030748402086, -1.0218124717262744, -0.8959221068158404, -0.7734940314247856, -0.6671434935220023, -0.5336857733899063, -0.3816574755357412, -0.29161087673834135, -0.17309384969561706, -0.05339460368627204, 0.06848024347088776, 0.15828211659462232, 0.2280793055415327, 0.35431338147484587, 0.43480510627034874, 0.4986527093096154, 0.5537304240173895, 0.6261102964754179, 0.6957256434779958, 0.7148716373433489, 0.7582185985675071, 0.8397185904494259, 0.8527817190580521, 0.8897458076340921, &
             1.3485557178362901, 1.4100564520861163, 1.5245130094695432, 1.5394854787019736, 1.5163778281947242, 1.3749694314333056, 1.2160585849035228, 0.9882317246213438, 0.7287282748280057, 0.4456317909919142, 0.174329001883184, -0.08314249078062917, -0.3703318878834143, -0.5197731810264961, -0.7456307038228134, -0.8983163369775448, -1.0210609293537283, -1.1325694312261778, -1.1872167191115837, -1.2303128662565612, -1.1635752906608618, -1.1366181255439949, -1.1143325488507188, -0.9566453631775117, -0.8555932167454776, -0.7255287019236673, -0.6054863016081078, -0.43728535957989434, -0.2276289420360999, -0.1244663569630013, 0.05384419703764767, 0.20686011064382656, 0.4031365843389501, 0.5579847371235589, 0.7033174462017056, 0.9474256754224659, 1.1190162927182088, 1.2326029774857028, 1.4189049778315361, 1.5297695457498912, &
             -1.9758514216280674, -1.6811036356993196, -1.4393679693100887, -1.0779006082963571, -0.7337931912174669, -0.2189941753254483, 0.1919445708944903, 0.5944234470589098, 0.850393105516085, 1.049798008613823, 1.2142668711181763, 1.2033705614127117, 1.230059224827435, 1.0210177131752327, 0.8880614290886265, 0.6940040971521174, 0.41674780998241634, 0.28161995443668403, 0.03424891416129952, -0.21795586614429763, -0.49862742114563025, -0.6619305026509881, -0.7207443580780598, -0.958937614573929, -1.0475571352581883, -1.0416931141536423, -1.0324618281900222, -1.0234784580976946, -1.0242682204528888, -0.8461579853478617, -0.6722214911029909, -0.492966559623527, -0.321261126282224, -0.05861948213877818, 0.27766540139242274, 0.41732701183840876, 0.7162331711558283, 1.1265674038691615, 1.4146082539850684, 1.804040045862745, &
             -2.1407006785624048, -1.0805568735683906, -0.3042327229719288, 0.6219138263143172, 1.014369869987673, 1.4247981784787198, 1.4744565229773687, 1.3330145845685064, 1.0266801133180212, 0.6621227879400782, 0.4072145153536709, 0.048448919497313864, -0.347740341476229, -0.613965371154217, -0.7932421998441845, -0.9129898751662757, -1.001562045972788, -0.8950071914103569, -0.7502386374091963, -0.6047307809868577, -0.5035176671798316, -0.25767631875524294, -0.04335890882119827, 0.13210094189985852, 0.39805291096418743, 0.46814712441809303, 0.735168968887883, 0.8279433128235392, 0.8522288866099015, 1.0309133422097239, 0.9689458683202629, 0.8194376351381301, 0.6292577529504737, 0.4942218318697948, 0.27043352358098177, -0.07164332919512675, -0.4515718482489549, -0.6908059434527288, -1.1764388578977747, -1.7174398416904961, &
             0.7797317676726162, -0.08664457661815766, 0.17762108863183143, -1.6850168891268635, -0.2488936349632652, -1.58154085759271, 0.13019531586723776, 0.3413195742419903, 0.694807351936538, 0.764142606603787, -0.22099329691615158, -0.08489271997893956, 1.373831668412292, -0.33547208236996573, 0.4889000587706866, 0.22324057674033926, 1.5063675253995814, 0.958968367113776, -1.8047184776647012, -1.9987438304179845, -0.37864689222689274, -0.546031276441644, -0.5693509852962297, -0.40477414545342305, -0.023584333562918964, 0.24150730223232975, -0.3620826334559182, 0.2569925172662072, -0.358807141373, 0.5763459916929877, -0.7263886505672454, 0.7279062932191324, 1.9322848706008275, 1.6192918141669006, -1.119038630610549, 1.7043430404044164, 1.085295356565518, -0.5939155123463035, -1.7671705066387882, -1.2550785150485222, &
             -0.8901818989508202, 0.6973654265739638, 2.4753961213337554, 1.0002855921359508, 1.720372003412172, -0.04258622882703431, -0.0044795598422154945, -0.4781913257056227, -0.7755737484188217, -0.749694423912793, -1.5906103874609656, -1.3831577053090143, 0.4086942926908348, -0.915606239943065, 0.13018087744810214, 0.34936926030053583, 1.2670317815742542, 1.1764959593387958, -0.05929056359246125, 0.08695674596625952, 0.5426310286437914, 0.5221277607921038, 1.0768122925042434, 0.0678286071017401, -0.06720639870241321, 0.6625768568134154, -0.2354407017105176, -0.44542978136637346, -1.3918143023152778, -0.6849575091709135, -1.5635037532196359, -0.1716428384041839, 0.2936786929427412, 0.12025637614547743, -0.7428196487835395, 0.3684259857534394, 0.2909633437999136, 0.15185730487009844, -0.2704534679820673, 1.4077039616576204, &
             0.425300769706358, -2.2067306209961535, -0.6985567915940346, -0.2523492550292019, 1.856231598451146, 0.06680903334829205, 0.7317872559635507, 0.2741527868166712, 0.2870720157513876, -0.11950196665836478, -0.5394222240857013, -0.4768856753534757, -1.7730239892459483, 0.990132270289617, -0.14517981090614943, -0.16811178925150427, 0.6952132674569004, -0.9388741511201266, 0.18945603201243078, -0.2816405145590186, 0.599160660601693, 0.9078289265561135, -0.8220550319996731, 0.41878398351194995, 0.8739735744175038, 0.5628782427648087, 0.13498480414564104, -0.1887128146937388, 1.232823279434726, 0.31452639777971625, -0.6826179726628224, -1.8601146271033715, 0.040429712783082716, -1.073055284773785, -2.9668789231327772, 0.5349921490649513, 2.0376287971602527, 0.08478924860009117, 0.8969419736038442, 0.07994172956890382, &
             0.6884939174387944, 0.7611658691201613, 0.8671621138068759, 0.8976311283816959, -0.03810995247223101, 0.3950341754447486, -1.1036063428657792, -0.8689588667156631, -1.240920893534697, -0.7139509019819348, -0.3116361201006657, -0.3178819356787596, -0.35745875002059974, 2.254607330516259, 0.9245446953402023, 1.0427771392916694, 0.9784777450612434, -0.11607772785846823, 0.7307771194089876, 0.2964625428668628, -0.11610293325694838, -0.3799926050420414, -2.340231770277648, -1.1213058479169322, -0.33207773884496145, -1.6996857377903072, -0.5487728714818453, -0.5971769913222891, 1.3508336835717878, 0.6072825923939345, 1.089557246757408, 0.35730653849216726, 0.7343379038308203, -0.08717943530437702, 0.40033418258495496, 0.6466128051046623, 0.48487465670495034, 0.4904361338831458, 0.17363072071820826, -2.2589459410292876, &
             -0.1124223621749161, -0.44051408218772437, -2.1217867972327578, 0.31261795818321997, 0.28996817271622044, 0.583780159804289, 1.9158402052509966, 1.4434215743346879, 0.05595770549451445, -0.9552354728153566, -1.1429509207147364, -1.0806799064495245, -0.8379535089622175, -1.1927267953337148, -0.607364248993462, -0.24502205481910186, 0.5986476612486711, 1.9750277619863885, 0.7572784336992865, 1.0480201864039516, 1.3034950375996952, 0.7556508203526529, -0.7124070872738679, 0.012943219005871226, 0.4911939761041664, -1.9274384661092128, -0.9502159144246344, -0.8118357269516286, -0.8093686725586665, -0.6669694630164572, -0.7602744610651714, 1.0462867196618815, 0.006687789693338572, 0.7638780537878009, 0.8903386202048241, 1.2134171858258596, -0.6462913542684801, 0.6282646495653872, 0.28407002895737754, -1.2700296498354382]
    contains
    procedure :: Init => TPCAReionization_Init
    procedure :: x_e => TPCAReionization_xe
    procedure :: get_timesteps => TPCAReionization_get_timesteps
    procedure :: ReadParams => TPCAReionization_ReadParams
    procedure, nopass :: SelfPointer => TPCAReionization_SelfPointer
    end type TPCAReionization

    public TBaseTauWithHeReionization, TTanhReionization, TExpReionization, TPCAReionization
    contains

    subroutine TBaseTauWithHeReionization_Init(this, State)
    use constants
    use MathUtils
    class(TBaseTauWithHeReionization) :: this
    class(TCAMBdata), target :: State
    procedure(obj_function) :: dtauda

    select type (State)
    class is (CAMBdata)
        this%State => State

        this%fHe =  State%CP%YHe/(mass_ratio_He_H*(1.d0-State%CP%YHe))
        if (this%Reionization) then

            if (this%optical_depth /= 0._dl .and. .not. this%use_optical_depth) &
                write (*,*) 'WARNING: You seem to have set the optical depth, but use_optical_depth = F'

            if (this%use_optical_depth.and.this%optical_depth<0.001 &
                .or. .not.this%use_optical_depth .and. this%Redshift<0.001) then
                this%Reionization = .false.
            end if

        end if

        if (this%Reionization) then

            if (this%fraction==Reionization_DefFraction) &
                this%fraction = 1._dl + this%fHe  !H + singly ionized He

            if (this%use_optical_depth) then
                call this%zreFromOptDepth()
                if (global_error_flag/=0) return
                if (FeedbackLevel > 0) write(*,'("Reion redshift       =  ",f6.3)') this%redshift
            end if

            call this%SetParamsForZre()

            !this is a check, agrees very well in default parameterization
            if (FeedbackLevel > 1) write(*,'("Integrated opt depth = ",f7.4)') this%State%GetReionizationOptDepth()

        end if
    end select
    end subroutine TBaseTauWithHeReionization_Init

    function TBaseTauWithHeReionization_SecondHelium_xe(this, z) result(xe)
    class(TBaseTauWithHeReionization) :: this
    real(dl), intent(in) :: z
    real(dl) xe, tgh, xod

    if (this%include_helium_fullreion .and. z < this%helium_redshiftstart) then
        !Effect of Helium becoming fully ionized is small so details not important
        xod = (this%helium_redshift - z)/this%helium_delta_redshift
        if (xod > 100) then
            tgh=1.d0
        else
            tgh=tanh(xod)
        end if

        xe = this%fHe*(tgh+1._dl)/2._dl
    else
        xe = 0.d0
    end if

    end function TBaseTauWithHeReionization_SecondHelium_xe


    subroutine TBaseTauWithHeReionization_ReadParams(this, Ini)
    use IniObjects
    class(TBaseTauWithHeReionization) :: this
    class(TIniFile), intent(in) :: Ini

    this%Reionization = Ini%Read_Logical('reionization')
    if (this%Reionization) then

        this%use_optical_depth = Ini%Read_Logical('re_use_optical_depth')

        if (this%use_optical_depth) then
            this%optical_depth = Ini%Read_Double('re_optical_depth')
        else
            this%redshift = Ini%Read_Double('re_redshift')
        end if

        call Ini%Read('re_ionization_frac',this%fraction)
        call Ini%Read('re_helium_redshift',this%helium_redshift)
        call Ini%Read('re_helium_delta_redshift',this%helium_delta_redshift)

        this%helium_redshiftstart  = Ini%Read_Double('re_helium_redshiftstart', &
            this%helium_redshift + 5*this%helium_delta_redshift)

    end if

    end subroutine TBaseTauWithHeReionization_ReadParams


    subroutine TBaseTauWithHeReionization_SetParamsForZre(this)
    class(TBaseTauWithHeReionization) :: this

    end subroutine TBaseTauWithHeReionization_SetParamsForZre

    subroutine TBaseTauWithHeReionization_Validate(this, OK)
    class(TBaseTauWithHeReionization),intent(in) :: this
    logical, intent(inout) :: OK

    if (this%Reionization) then
        if (this%use_optical_depth) then
            if (this%optical_depth<0 .or. this%optical_depth > 0.9  .or. &
                this%include_helium_fullreion .and. this%optical_depth<0.01) then
                OK = .false.
                write(*,*) 'Optical depth is strange. You have:', this%optical_depth
            end if
        end if
        if (this%fraction/= Reionization_DefFraction .and. (this%fraction < 0 .or. this%fraction > 1.5)) then
            OK = .false.
            write(*,*) 'Reionization fraction strange. You have: ',this%fraction
        end if
    end if

    end subroutine TBaseTauWithHeReionization_Validate

    subroutine TBaseTauWithHeReionization_zreFromOptDepth(this)
    !General routine to find zre parameter given optical depth
    class(TBaseTauWithHeReionization) :: this
    real(dl) try_b, try_t
    real(dl) tau, last_top, last_bot
    integer i

    try_b = this%min_redshift
    try_t = this%max_redshift
    i=0
    do
        i=i+1
        this%redshift = (try_t + try_b)/2
        call this%SetParamsForZre()
        tau = this%State%GetReionizationOptDepth()

        if (tau > this%optical_depth) then
            try_t = this%redshift
            last_top = tau
        else
            try_b = this%redshift
            last_bot = tau
        end if
        if (abs(try_b - try_t) < 1e-2_dl/this%tau_solve_accuracy_boost) then
            if (try_b==this%min_redshift) last_bot = this%min_redshift
            if (try_t/=this%max_redshift) this%redshift  = &
                (try_t*(this%optical_depth-last_bot) + try_b*(last_top-this%optical_depth))/(last_top-last_bot)
            exit
        end if
        if (i>100) call GlobalError('TBaseTauWithHeReionization_zreFromOptDepth: failed to converge',error_reionization)
    end do

    if (abs(tau - this%optical_depth) > 0.002 .and. global_error_flag==0) then
        write (*,*) 'TBaseTauWithHeReionization_zreFromOptDepth: Did not converge to optical depth'
        write (*,*) 'tau =',tau, 'optical_depth = ', this%optical_depth
        write (*,*) try_t, try_b
        write (*,*) '(If running a chain, have you put a constraint on tau?)'
        call GlobalError('Reionization did not converge to optical depth',error_reionization)
    end if

    end subroutine TBaseTauWithHeReionization_zreFromOptDepth

    real(dl) function TBaseTauWithHeReionization_GetZreFromTau(P, tau)
    type(CAMBparams) :: P, P2
    real(dl) tau
    integer error
    type(CAMBdata) :: State

    P2 = P

    select type(Reion=>P2%Reion)
    class is (TBaseTauWithHeReionization)
        Reion%Reionization = .true.
        Reion%use_optical_depth = .true.
        Reion%optical_depth = tau
    end select
    call State%SetParams(P2,error)
    if (error/=0)  then
        TBaseTauWithHeReionization_GetZreFromTau = -1
    else
        select type(Reion=>State%CP%Reion)
        class is (TBaseTauWithHeReionization)
            TBaseTauWithHeReionization_GetZreFromTau = Reion%redshift
        end select
    end if

    end function  TBaseTauWithHeReionization_GetZreFromTau

    function TTanhReionization_xe(this, z, tau, xe_recomb)
    !a and time tau are redundant, both provided for convenience
    !xe_recomb is xe(tau_start) from recombination (typically very small, ~2e-4)
    !xe should map smoothly onto xe_recomb
    class(TTanhReionization) :: this
    real(dl), intent(in) :: z
    real(dl), intent(in), optional :: tau, xe_recomb
    real(dl) TTanhReionization_xe
    real(dl) tgh, xod
    real(dl) xstart

    xstart = PresentDefault(0._dl, xe_recomb)

    xod = (this%WindowVarMid - (1+z)**Tanh_zexp)/this%WindowVarDelta
    if (xod > 100) then
        tgh=1.d0
    else
        tgh=tanh(xod)
    end if
    TTanhReionization_xe =(this%fraction-xstart)*(tgh+1._dl)/2._dl+xstart + &
        this%SecondHelium_xe(z)

    end function TTanhReionization_xe

    subroutine TTanhReionization_get_timesteps(this, n_steps, z_start, z_complete)
    !minimum number of time steps to use between tau_start and tau_complete
    !Scaled by AccuracyBoost later
    !steps may be set smaller than this anyway
    class(TTanhReionization) :: this
    integer, intent(out) :: n_steps
    real(dl), intent(out):: z_start, z_Complete

    n_steps = nint(50 * this%timestep_boost)
    z_start = this%redshift + this%delta_redshift*8
    z_complete = max(0.d0,this%redshift-this%delta_redshift*8)

    end subroutine TTanhReionization_get_timesteps

    subroutine TTanhReionization_SetParamsForZre(this)
    class(TTanhReionization) :: this

    this%WindowVarMid = (1._dl+this%redshift)**Tanh_zexp
    this%WindowVarDelta = Tanh_zexp*(1._dl+this%redshift)**(Tanh_zexp-1._dl)*this%delta_redshift

    end subroutine TTanhReionization_SetParamsForZre

    subroutine TTanhReionization_ReadParams(this, Ini)
    use IniObjects
    class(TTanhReionization) :: this
    class(TIniFile), intent(in) :: Ini

    call this%TBaseTauWithHeReionization%ReadParams(Ini)
    if (this%Reionization) call Ini%Read('re_delta_redshift',this%delta_redshift)

    end subroutine TTanhReionization_ReadParams

    subroutine TTanhReionization_Validate(this, OK)
    class(TTanhReionization),intent(in) :: this
    logical, intent(inout) :: OK

    call this%TBaseTauWithHeReionization%Validate(OK)
    if (this%Reionization) then
        if (.not. this%use_optical_depth) then
            if (this%redshift < 0 .or. this%Redshift +this%delta_redshift*3 > this%max_redshift .or. &
                this%include_helium_fullreion .and. this%redshift < this%helium_redshift) then
                OK = .false.
                write(*,*) 'Reionization redshift strange. You have: ',this%Redshift
            end if
        end if
        if (this%delta_redshift > 3 .or. this%delta_redshift<0.1 ) then
            !Very narrow windows likely to cause problems in interpolation etc.
            !Very broad likely to conflict with quasar data at z=6
            OK = .false.
            write(*,*) 'Reionization delta_redshift is strange. You have: ',this%delta_redshift
        end if
    end if

    end subroutine TTanhReionization_Validate

    subroutine TTanhReionization_SelfPointer(cptr,P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type (TTanhReionization), pointer :: PType
    class (TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType

    end subroutine TTanhReionization_SelfPointer

    subroutine TExpReionization_Init(this, State)
    class(TExpReionization) :: this
    class(TCAMBdata), target :: State

    this%min_redshift = this%reion_redshift_complete
    call this%TBaseTauWithHeReionization%Init(State)

    end subroutine TExpReionization_Init


    function TExpReionization_xe(this, z, tau, xe_recomb)
    !a and time tau are redundant, both provided for convenience
    !xe_recomb is xe(tau_start) from recombination (typically very small, ~2e-4)
    !xe should map smoothly onto xe_recomb
    class(TExpReionization) :: this
    real(dl), intent(in) :: z
    real(dl), intent(in), optional :: tau, xe_recomb
    real(dl) TExpReionization_xe
    real(dl) lam, xstart, smoothing

    xstart = PresentDefault(0._dl, xe_recomb)

    if (z <= this%reion_redshift_complete + 1d-6) then
        TExpReionization_xe = this%fraction
    else
        lam = -log(0.5)/(this%redshift - this%reion_redshift_complete)**this%reion_exp_power
        smoothing = 1/(1+this%reion_exp_smooth_width/(z-this%reion_redshift_complete)**2)
        TExpReionization_xe = exp(-lam*(z-this%reion_redshift_complete)**this%reion_exp_power*smoothing) &
            *(this%fraction-xstart) + xstart
    end if

    TExpReionization_xe = TExpReionization_xe +  this%SecondHelium_xe(z)

    end function TExpReionization_xe

    subroutine TExpReionization_get_timesteps(this, n_steps, z_start, z_complete)
    !minimum number of time steps to use between tau_start and tau_complete
    !Scaled by AccuracyBoost later
    !steps may be set smaller than this anyway
    class(TExpReionization) :: this
    integer, intent(out) :: n_steps
    real(dl), intent(out):: z_start, z_complete
    real(dl) lam

    n_steps = nint(50 * this%timestep_boost)
    lam = -log(0.5)/(this%redshift - this%reion_redshift_complete)**this%reion_exp_power
    z_start = this%reion_redshift_complete  + (-log(0.0001)/lam)**(1/this%reion_exp_power)
    z_complete = this%reion_redshift_complete

    end subroutine TExpReionization_get_timesteps

    subroutine TExpReionization_ReadParams(this, Ini)
    use IniObjects
    class(TExpReionization) :: this
    class(TIniFile), intent(in) :: Ini

    call this%TBaseTauWithHeReionization%ReadParams(Ini)
    if (this%Reionization)then
        call Ini%Read('reion_redshift_complete',this%reion_redshift_complete)
        call Ini%Read('reion_exp_smooth_width',this%reion_exp_smooth_width)
        call Ini%Read('reion_exp_power',this%reion_exp_power)
    end if

    end subroutine TExpReionization_ReadParams

    subroutine TExpReionization_SelfPointer(cptr,P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type (TExpReionization), pointer :: PType
    class (TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType

    end subroutine TExpReionization_SelfPointer



    !MT PCA reionization
    subroutine TPCAReionization_Init(this, State)
    class(TPCAReionization) :: this
    class(TCAMBdata), target :: State
    integer i

    call this%TBaseTauWithHeReionization%Init(State)

    end subroutine TPCAReionization_Init

    function TPCAReionization_xe(this, z, tau, xe_recomb)
    !a and time tau are redundant, both provided for convenience
    !xe_recomb is xe(tau_start) from recombination (typically very small, ~2e-4)
    !xe should map smoothly onto xe_recomb
    class(TPCAReionization) :: this
    real(dl), intent(in) :: z
    real(dl), intent(in), optional :: tau, xe_recomb
    real(dl) TPCAReionization_xe
    real(dl) interp_scale, xstart
    real(dl) zstart, zend
    real, dimension(10,40) :: pca
    real, dimension(10) :: modes
    integer iz, imode, n

    xstart = PresentDefault(0._dl, xe_recomb)
!    call this%get_timesteps(n, zstart, zend)
    zstart = MAXVAL( this%reion_redshift)
    zend   = MINVAL( this%reion_redshift)

    pca = reshape( this%reion_pca, shape(pca), order=(/2,1/))
    modes = [ &
         this%reion_pca_mode1, this%reion_pca_mode2, this%reion_pca_mode3, this%reion_pca_mode4, this%reion_pca_mode5, &
         this%reion_pca_mode6, this%reion_pca_mode7, this%reion_pca_mode8, this%reion_pca_mode9, this%reion_pca_mode10]

    TPCAReionization_xe = xstart
    if (z < zstart .and. z > zend) then
       TPCAReionization_xe = TPCAReionization_xe + 0.15

       iz = 1
       do while( this%reion_redshift(iz) < z)
          iz = iz+1
       end do
       interp_scale = (z - this%reion_redshift(iz-1)) / (this%reion_redshift(iz) - this%reion_redshift(iz-1))

       do imode = 1, size(modes)
          TPCAReionization_xe = TPCAReionization_xe + (this%fraction-xstart) * modes(imode) * (pca(imode,iz-1) + (pca(imode,iz) - pca(imode,iz-1))*interp_scale)
       end do

    else if (z < zend) then
       TPCAReionization_xe = this%fraction
    end if

    TPCAReionization_xe = TPCAReionization_xe + this%SecondHelium_xe(z)

    end function TPCAReionization_xe

    subroutine TPCAReionization_get_timesteps(this, n_steps, z_start, z_complete)
    !minimum number of time steps to use between tau_start and tau_complete
    !Scaled by AccuracyBoost later
    !steps may be set smaller than this anyway
    class(TPCAReionization) :: this
    integer, intent(out) :: n_steps
    real(dl), intent(out):: z_start, z_complete
    real(dl) lam

    n_steps = nint(50 * this%timestep_boost)
    z_start = 30
    z_complete = 5

    end subroutine TPCAReionization_get_timesteps

    subroutine TPCAReionization_ReadParams(this, Ini)
    use IniObjects
    class(TPCAReionization) :: this
    class(TIniFile), intent(in) :: Ini
    integer i

!    call this%TBaseTauWithHeReionization%ReadParams(Ini)
    this%Reionization = Ini%Read_Logical('reionization')
    if (this%Reionization) then

       this%use_optical_depth = Ini%Read_Logical('re_use_optical_depth')
       if (this%use_optical_depth) then
          call GlobalError('TPCAReionization_ReadParams: Cannot use optical_depth for PCA model',error_reionization)
       end if
       
       call Ini%Read('re_ionization_frac',this%fraction)
       call Ini%Read('re_helium_redshift',this%helium_redshift)
       call Ini%Read('re_helium_delta_redshift',this%helium_delta_redshift)
       
       this%helium_redshiftstart  = Ini%Read_Double('re_helium_redshiftstart', &
            this%helium_redshift + 5*this%helium_delta_redshift)
       
       call Ini%Read('reion_pca_mode1',this%reion_pca_mode1)
       call Ini%Read('reion_pca_mode2',this%reion_pca_mode2)
       call Ini%Read('reion_pca_mode3',this%reion_pca_mode3)
       call Ini%Read('reion_pca_mode4',this%reion_pca_mode4)
       call Ini%Read('reion_pca_mode5',this%reion_pca_mode5)
       call Ini%Read('reion_pca_mode6',this%reion_pca_mode6)
       call Ini%Read('reion_pca_mode7',this%reion_pca_mode7)
       call Ini%Read('reion_pca_mode8',this%reion_pca_mode8)
       call Ini%Read('reion_pca_mode9',this%reion_pca_mode9)
       call Ini%Read('reion_pca_mode10',this%reion_pca_mode10)
    end if

    end subroutine TPCAReionization_ReadParams

    subroutine TPCAReionization_SelfPointer(cptr,P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type (TPCAReionization), pointer :: PType
    class (TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType

    end subroutine TPCAReionization_SelfPointer




    end module Reionization
