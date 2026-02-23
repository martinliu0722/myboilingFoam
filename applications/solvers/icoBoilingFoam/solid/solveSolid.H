if (finalIter)
{
    mesh.data::add("finalIteration", true);
}

{
    // 查找之前准备好的源项场
    // 确保在 createSolidFields.H 中已经创建了这个场，
    // 并且在 icoBoilingFoam.C 主循环的映射步骤中已经更新了它。
    const volScalarField& microSource = 
        mesh.lookupObject<volScalarField>("solidMicrolayerSource");
    
    for (int nonOrth = 0; nonOrth <= nNonOrthCorr; ++nonOrth)
    {
        fvScalarMatrix hEqn
        (
            fvm::ddt(betav*rho, h)
          - (
               thermo.isotropic()
             ? fvm::laplacian(betav*thermo.alpha(), h, "laplacian(alpha,h)")
             : fvm::laplacian(betav*taniAlpha(), h, "laplacian(alpha,h)")
            )
          ==
            fvOptions(rho, h)
            // 添加微液层吸热源项
            // 注意符号：
            // 如果 microSource 在映射时定义为负值 (Sink, W/m3)，则用 + 号
            // 如果 microSource 定义为正值 (代表流体带走的热量)，则用 - 号
            + microSource
        );

        hEqn.relax();

        fvOptions.constrain(hEqn);

        hEqn.solve(mesh.solver(h.select(finalIter)));

        fvOptions.correct(h);
    }

    thermo.correct();

    Info<< "Min/max T:" << min(thermo.T()).value() << ' '
        << max(thermo.T()).value() << endl;
}

if(computeFluxes)
{
    #include "computeSolidHeatFluxes.H"
}

if (finalIter)
{
    mesh.data::remove("finalIteration");
}
