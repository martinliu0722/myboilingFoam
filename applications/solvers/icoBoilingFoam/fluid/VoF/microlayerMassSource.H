//- Update mDot to account for microlayer models at the walls

if(mixture.found("microlayerModel"))
{
    //Info<< "\nStarting microlayer stuff"<< endl;
    dictionary microDict(mixture.subDict("microlayerModel"));
    Switch debugML(mixture.lookupOrDefault<Switch>("debug",false));

    //- Coefficients for polynomium
    scalarField coeffs(microDict.lookup("coefficients"));

    if(coeffs.size() != 7)
    {
        Info<<"\nThe microlayer model requires 7 coefficients!\n"
            << abort(FatalError);
    }

    wordList activePatches(microDict.lookup("patches"));

    scalarField ap(volPointInterpolation::New(mesh).interpolate(alpha1));

    cutFaceIso faceCut(mesh,ap);

    forAll(activePatches,patchi)
    {
        label patchID(mesh.boundaryMesh().findPatchID(activePatches[patchi]));

        if(patchID == -1)
        {
            continue;
        }

        const fvPatchScalarField& TBf(T.boundaryFieldRef()[patchID]);
        const fvPatchScalarField& aBf(alpha1.boundaryFieldRef()[patchID]);
        scalarField& kappaBf(kappa.boundaryFieldRef()[patchID]);
        const labelUList& fCells(TBf.patch().faceCells());
        const scalarField& deltas(TBf.patch().deltaCoeffs());
        const scalarField pA(mesh.magSf().boundaryField()[patchID]);

        //- Cache wall gradient
        scalarField dTdn(TBf.snGrad());

        const label startFace(TBf.patch().start());

        // --- [DEBUG START] initialize variables ---
        // label nContactLines = 0;           // 统计接触线个数
        // scalar minLc = GREAT;              // 初始化为极大值
        // scalar maxLc = -GREAT;             // 初始化为极小值
        // scalar sumLc = 0.0;                // 用于计算平均值
        // --- [DEBUG END] ---

        forAll(aBf,facei)
        {
            label cellI = fCells[facei];

            LatTSource[cellI] = 0;
            mDot0[cellI] = 0;

            //- Select interface cells only
            if(aBf[facei] < 1e-3 || aBf[facei] > 0.999)
            {
                continue;
            }

            if(faceCut.calcSubFace(startFace + facei,0.5) != 0)
            {
                continue;
            }

            if(faceCut.subFacePoints().size() == 0)
            {
                continue;
            }
            // Info<< "\npoints = "<< faceCut.surfacePoints() <<  endl;

            //- Compute the contact length
            point p0(faceCut.surfacePoints()[0]);
            point p1(faceCut.surfacePoints()[1]);

            // 1. 安全计算 lc (确保一定是标量数值)
            scalar lcValue = mag(p0 - p1);

            /*
            // --- [DEBUG START] 收集数据 ---
            nContactLines++;
            if (lcValue < minLc) minLc = lcValue;
            if (lcValue > maxLc) maxLc = lcValue;
            sumLc += lcValue;
            // --- [DEBUG END] ---
            */

            dimensionedScalar lc
            (
                "lc",
                dimLength,
                lcValue
            );

            //- Set conduction to zero.
            //  This does not affect cht, which is based on the true kappa.
            kappaBf[facei] = 0.;

            scalar deltaTmp(TBf[facei] - Tsat.value());

            //- Re-compute mDot0 at the interface
            mDot0[cellI] +=
            (
                lcValue/(hf.value()*mesh.V()[cellI])
                *
                (
                //     (coeffs[0])
                //   + (coeffs[1] * pow(deltaTmp, 1))
                //   + (coeffs[2] * pow(deltaTmp, 2))
                    
                //     (coeffs[0] * deltaTmp)
                //   + (coeffs[1] * pow(deltaTmp, 2))
                //   + (coeffs[2] * pow(deltaTmp, 3))
                //   + (coeffs[3] * pow(deltaTmp, 4))
                    
                       coeffs[0]
                    + (coeffs[1] * std::pow(deltaTmp, 1.0))
                    + (coeffs[2] * std::pow(deltaTmp, 2.0))
                    + (coeffs[3] * std::pow(deltaTmp, 3.0))
                    + (coeffs[4] * std::pow(deltaTmp, 4.0))
                    + (coeffs[5] * std::pow(deltaTmp, 5.0))
                    + (coeffs[6] * std::pow(deltaTmp, 6.0))

                    // (coeffs[0] * deltaTmp)
                    // + (coeffs[1] * std::pow(deltaTmp, 1.0 / 2.0))
                    // + (coeffs[2] * std::pow(deltaTmp, 1.0 / 3.0))
                    // + (coeffs[3] * std::pow(deltaTmp, 1.0 / 4.0))
                    // + (coeffs[4] * std::pow(deltaTmp, 1.0 / 5.0))
                    // + (coeffs[5] * std::pow(deltaTmp, 1.0 / 6.0))
                    // + (coeffs[6] * std::pow(deltaTmp, 1.0 / 7.0))
                    // + (coeffs[7] * std::pow(deltaTmp, 1.0 / 8.0))
                )
            );

            LatTSource[cellI] +=
            (
                // lcValue
                -lcValue/(mesh.V()[cellI])
                *
                (
                       coeffs[0]
                    + (coeffs[1] * std::pow(deltaTmp, 1.0))
                    + (coeffs[2] * std::pow(deltaTmp, 2.0))
                    + (coeffs[3] * std::pow(deltaTmp, 3.0))
                    + (coeffs[4] * std::pow(deltaTmp, 4.0))
                    + (coeffs[5] * std::pow(deltaTmp, 5.0))
                    + (coeffs[6] * std::pow(deltaTmp, 6.0))
                )
            );

            if(debugML)
            {
                fvPatchScalarField& lcBf
                (
                    lcFluid[i].boundaryFieldRef()[patchID]
                );

                lcBf[facei] = lc.value();
            }
        }

        /*
        // --- [DEBUG START] 并行归约与输出 ---
        reduce(nContactLines, sumOp<label>());
        reduce(minLc, minOp<scalar>());
        reduce(maxLc, maxOp<scalar>());
        reduce(sumLc, sumOp<scalar>());

        if (nContactLines > 0)
        {
            Info << "\n[DEBUG LC Stats] Patch: " << activePatches[patchi] << nl
                << "    Count   : " << nContactLines << nl
                << "    Min lc  : " << minLc << " m" << nl
                << "    Max lc  : " << maxLc << " m" << nl
                << "    Avg lc  : " << sumLc / nContactLines << " m" << nl
                << "    Total L : " << sumLc << " m" << endl;
        }
        else
        {
            Info << "\n[DEBUG LC Stats] Patch: " << activePatches[patchi] 
                << " -> No contact lines found (count = 0)." << endl;
        }
        // --- [DEBUG END] ---
        */

        // ========================================================================
        // 将microlayer的 LatTSource 映射到固体区域
        // ========================================================================
/*
                    const polyPatch &pp = mesh.boundaryMesh()[patchID];

        // 依然需要检查这是否是一个映射边界，否则无法传数据
        if (isA<mappedPatchBase>(pp))
        {
            // Info << "\n------Find mapped Patch Base!------" << endl;

            const mappedPatchBase &mpp = refCast<const mappedPatchBase>(pp);

            // 1. 获取邻居固体网格
            const word &nbrRegionName = mpp.sampleRegion();
            const fvMesh &solidMesh = mesh.time().lookupObject<fvMesh>(nbrRegionName);

            // Info << "The sampleRegion name (solid) is: " << nbrRegionName << "\n" << endl;

            // 2. 获取固体目标场 (需确保它在固体区域已注册)
            // 这里用指针的方式定义变量,保证了solidMicrolayerSource会同步energySrcSolid的数值变化
            volScalarField &energySrcSolid = const_cast<volScalarField &>
            (
                solidMesh.lookupObject<volScalarField>("solidMicrolayerSource")
            );

            volScalarField &solidRefineField = const_cast<volScalarField &>
            (
                solidMesh.lookupObject<volScalarField>("solidrefine")
                );

            const fvPatchScalarField &alpha1Bf = alpha1.boundaryField()[patchID];

            scalarField alphaFluid = alpha1Bf; // 拷贝一份数据

            // 3. 计算流体侧边界热通量 q [W/m2]
            // 利用当前上下文中的 mesh, LatTSource
            scalarField qFluid(pp.size(), 0.0);
            const labelUList &fluidFaceCells = pp.faceCells();
            const scalarField &fluidMagSf = pp.magFaceAreas();

            forAll(fluidFaceCells, faceI)
            {
                label cellI = fluidFaceCells[faceI];
                // 这里的 LatTSource 刚刚在上面的循环里被赋值了，现在直接拿来用
                qFluid[faceI] = LatTSource[cellI] * mesh.V()[cellI] / fluidMagSf[faceI];
            }

            // 4. 数据传输 (Fluid Patch -> Solid Patch)
            scalarField qOnSolidPatch = qFluid;
            mpp.reverseDistribute(qOnSolidPatch);

            mpp.reverseDistribute(alphaFluid);

            
            // --- [DEBUG 代码开始] ---
            gMin/gMax 支持并行，能看到所有处理器的极值
            Info << "\n[DEBUG] qOnSolidPatch Stats:" << nl
                << "    Size: " << qOnSolidPatch.size() << nl
                << "    Min : " << gMin(qOnSolidPatch) << nl
                << "    Max : " << gMax(qOnSolidPatch) << nl
                << "    Avg : " << gAverage(qOnSolidPatch) << endl;
            // --- [DEBUG 代码结束] ---
            

            // 5. 固体侧赋值 (Face -> Owner Cell)
            const word &solidPatchID = mpp.samplePatch();
            const polyPatch &solidPatch = solidMesh.boundaryMesh()[solidPatchID];
            const scalarField &solidMagSf = solidPatch.magFaceAreas();
            const labelUList &solidBoundaryCells = solidPatch.faceCells();
            // Info << "The solidPatchID (solid) is: " << solidPatchID << "\n" << endl;

            forAll(solidBoundaryCells, faceI)
            {
                label cellID = solidBoundaryCells[faceI];
                scalar cellVol = solidMesh.V()[cellID];
                scalar faceArea = solidMagSf[faceI];
                scalar mappedAlpha = alphaFluid[faceI];

                // 流体放热 = 固体吸热
                scalar flux = qOnSolidPatch[faceI];

                // 累加源项 [W/m3] 
                // 注意：这里用 += 还是 = 取决于你是否还有其他机制在写这个变量
                // 如果这是唯一的源，用 =；如果是累加，用 +=
                energySrcSolid[cellID] = flux * faceArea / cellVol;

                solidRefineField[cellID] = max(solidRefineField[cellID], mappedAlpha);
            }

            // 6. 更新并行边界(可选,没有测试是否有用)
            energySrcSolid.correctBoundaryConditions();
            solidRefineField.correctBoundaryConditions();

            // Debug 输出 (可选)
            // Info << "Mapped microlayer source from " << pp.name() 
            //      << " to " << nbrRegionName << "\n" << endl;
        

        }
*/
    }
}
