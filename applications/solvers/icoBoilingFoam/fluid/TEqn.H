{

    kappa = alpha1*k1 + alpha2*k2;
    surfaceScalarField kappaf(fvc::interpolate(kappa));

    // Generate mixture phase Cp
    volScalarField Cpmix
    (
        IOobject
        (
            "Cpmix",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alpha1*cp1 + alpha2*cp2,
        alpha1.boundaryField().types()
    );
    Cpmix.oldTime();

    // The sensible heat source MUST BE calculated before "computeMassSource.H"
    // or the wrong T field will be obtaied!!!
    // Assumption: vapor appear in Tsat
    volScalarField impSenTSource(-continuitySource * Cpmix);
    volScalarField expSenTSource(continuitySource *Cpmix * Tsat);

    #include "computeMassSource.H"

    volScalarField impTSource(impLatTSource + impSenTSource);
    volScalarField expTSource(expLatTSource + expSenTSource);

    fvScalarMatrix TEqn
     (
         fvm::ddt(rhoCp,T)
       + fvm::div(rhoCpPhi,T)
       - fvm::Sp(fvc::ddt(rhoCp) + fvc::div(rhoCpPhi), T)
       - fvm::laplacian(kappaf,T)
       ==
         fvm::Sp(impTSource, T)
       + expTSource
     );

     TEqn.relax();
     TEqn.solve();

     Info<< "min/max(T) = " << min(T).value() << ", " << max(T).value() <<endl;

     //update continuity source with new T
    #include "computeMassSource.H"

    // ---------------------------------------------------------
    // 计算流体侧壁面热通量 qwFluid
    // ---------------------------------------------------------

    if (mixture.found("writeQwall"))
    {
        // 1. 获取场引用
        volVectorField& qw = qwFluid[i];
        const volScalarField& LatT = LatTSourceFluid[i];
        const volScalarField& kappa = kappaFluid[i];
        const volScalarField& T = TFluid[i];
        
        // 2. 获取网格引用
        const fvMesh& mesh = fluidRegions[i];
        
        // 获取网格中心坐标场
        const volVectorField& C = mesh.C();

        // 3. 全局初始化：使用傅里叶导热定律
        qw = -kappa * fvc::grad(T);

        // 4. 边界修正循环
        const fvBoundaryMesh& patches = mesh.boundary();

        forAll(patches, patchI)
        {
            const fvPatch& patch = patches[patchI];
            
            // --- 核心修改：使用字符串匹配代替 isA<Type> ---
            // 这种方法不需要 include 任何头文件，编译绝对通过！
            // 只要满足以下任一条件，就认为是我们要处理的壁面：
            // 1. 类型是 wall
            // 2. 类型是 mappedWall (CHT 耦合面)
            // 3. 名字里包含 "fluid_to_solid" (你的特定需求)
            if (
                patch.type() == "wall" || 
                patch.type() == "mappedWall" ||
                patch.name() == "fluid_to_solid"
               )
            {
                // 获取该 Patch 的面法向量 (单位向量，指向流体外部)
                const vectorField& n = patch.nf();
                
                // 获取该 Patch 的面中心坐标
                const vectorField& Cf = patch.Cf();

                // 获取该 Patch 对应内部网格的索引列表
                const labelList& faceCells = patch.faceCells();

                // 遍历该 Patch 上的所有面
                forAll(faceCells, faceI)
                {
                    // 获取对应的内部网格 ID
                    label cellI = faceCells[faceI];

                    // 核心判定：只有当该网格有潜热源项 (LatTSource != 0) 时才修正
                    if (mag(LatT[cellI]) > SMALL)
                    {
                        // A. 计算精确的 dz (几何投影法)
                        // 向量 d = 面中心(Cf) - 网格中心(C)
                        vector d = Cf[faceI] - C[cellI];
                        
                        // 将 d 投影到法向量 n 上，取绝对值
                        // 这是最准确的 "Mesh.dz"
                        scalar dz = 2.0*mag(d & n[faceI]);

                        // B. 计算热流密度大小 Q [W/m2] = Source [W/m3] * dz [m]
                        scalar qMag = LatT[cellI] * dz;

                        // C. 确定方向
                        // patch.nf() 指向外部(固体)。
                        // 热流从固体流向流体(向内)。
                        // 所以方向是 -n
                        vector dir = -n[faceI];

                        // D. 强制覆盖赋值
                        qw[cellI] = dir * qMag;
                    }
                }
            }
        }
        
        // 5. 修正边界场
        qw.correctBoundaryConditions();
        
        Info<< "    Corrected qwFluid direction on fluid_to_solid interface." << endl;
    }
}
