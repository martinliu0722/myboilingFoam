/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 Federico Municchi. All rights reserved.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "bubbleGenerator.H"
#include "Time.H"
#include "fvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "Pstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(bubbleGenerator, 0);
    addToRunTimeSelectionTable(functionObject, bubbleGenerator, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::bubbleGenerator::bubbleGenerator
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    runTime_(runTime),
    fixDT_(dict.lookupOrDefault<bool>("fixTimeStepAfterNucleation", false)),
    maxLifeSteps_(dict.lookupOrDefault<label>("maxLifeSteps", 1)),
    TName_("T"),
    alphaName_("alpha.liquid")
{
    read(dict);

    location_.reset
    (
        locationModel::New
        (
            dict,
            mesh_
        )
    );

    activation_.reset
    (
        activationModel::New
        (
            dict,
            mesh_
        )
    );

    nucleation_.reset
    (
        nucleationModel::New
        (
            dict,
            mesh_
        )
    );

    if(fixDT_)
    {
        newDT_ = readScalar(dict.lookup("deltaT"));
    }

    // [User Added] 读取自定义参数
    dict.readIfPresent("TemperatureField", TName_);
    dict.readIfPresent("AlphaField", alphaName_);
    dict.readIfPresent("maxLifeSteps", maxLifeSteps_);

    // [User Added] 读取成核点列表
    if (dict.found("nucleationSites"))
    {
        const PtrList<dictionary> siteDicts(dict.lookup("nucleationSites"));
        mySites_.setSize(siteDicts.size());

        forAll(siteDicts, i)
        {
            const dictionary& sd = siteDicts[i];
            sd.lookup("location") >> mySites_[i].location;
            mySites_[i].Tactivate = readScalar(sd.lookup("Tactivate"));
            mySites_[i].radius    = readScalar(sd.lookup("radius"));
            mySites_[i].waitTime  = sd.lookupOrDefault<scalar>("waitTime", 0.05);
            
            Info<< "  Added site " << i << ": pos=" << mySites_[i].location 
                << ", T_act=" << mySites_[i].Tactivate 
                << ", R=" << mySites_[i].radius << endl;
        }
    }

    readState();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::bubbleGenerator::read(const dictionary& dict)
{

    return true;
}



bool Foam::functionObjects::bubbleGenerator::execute()
{
    /* ---- Martin's codes come here ---- */
    /* ---- Version 2.0: Active sites with face temperature*/
    const fvMesh &mesh = mesh_;
    const volVectorField& C = mesh.C();
    const scalar currentTime = runTime_.value();
    
    // =========================================================
    // 1. 读取饱和温度 Tsat (兼容 constant/fluid 路径)
    // =========================================================
    scalar Tsat = 373.15; // 默认值
    
    // 显式指定在 constant/fluid 下查找
    IOobject transPropIO(
        "transportProperties",
        mesh.time().constant(), 
        "fluid", 
        mesh, 
        IOobject::READ_IF_PRESENT, 
        IOobject::NO_WRITE
    );

    if (transPropIO.typeHeaderOk<IOdictionary>(true))
    {
        IOdictionary transDict(transPropIO);
        if (transDict.found("PhaseChangeProperties")) 
        {
            transDict.subDict("PhaseChangeProperties").readIfPresent("Tsat", Tsat);
        }
    }
    else if (mesh.foundObject<IOdictionary>("transportProperties")) 
    {
        // 回退机制：在默认注册表中查找
        const IOdictionary& transDict = mesh.lookupObject<IOdictionary>("transportProperties");
        if (transDict.found("PhaseChangeProperties")) 
        {
            transDict.subDict("PhaseChangeProperties").readIfPresent("Tsat", Tsat);
        }
    }
    
    // =========================================================
    // 2. 获取 alpha 和 T 场
    // =========================================================
    volScalarField* ptrAlpha = const_cast<volScalarField*>(
        mesh.findObject<volScalarField>(alphaName_)
    );
    volScalarField* ptrT = const_cast<volScalarField*>(
        mesh.findObject<volScalarField>(TName_)
    );

    // 如果找不到场，安静地返回 false (或仅在 Master 警告一次，这里保持安静)
    if (!ptrAlpha || !ptrT) return false;

    volScalarField& alpha = *ptrAlpha;
    volScalarField& T = *ptrT;

    // =========================================================
    // 3. 准备流固耦合边界数据
    // =========================================================
    label patchID = mesh.boundaryMesh().findPatchID("fluid_to_solid");
    
    vectorField faceCenters; 
    const fvPatchScalarField* TpatchPtr = nullptr;

    if (patchID != -1)
    {
        const polyPatch& pp = mesh.boundaryMesh()[patchID];
        // 直接拷贝面心坐标，避免引用临时对象的编译错误
        faceCenters = pp.faceCentres(); 
        TpatchPtr = &T.boundaryField()[patchID];
    }

    // =========================================================
    // 4. 遍历成核点 (Local Search)
    // =========================================================
    List<bool> localTriggers(mySites_.size(), false);
    List<point> localBubbleCenters(mySites_.size(), vector::zero);

    forAll(mySites_, i)
    {
        NucleationSite& site = mySites_[i];
        
        if (!site.isActive)
        {
            scalar minDist = GREAT;
            label minFaceID = -1;

            // 本地搜索最近的面
            if (faceCenters.size() > 0)
            {
                forAll(faceCenters, facei)
                {
                    scalar d = mag(faceCenters[facei] - site.location);
                    if (d < minDist)
                    {
                        minDist = d;
                        minFaceID = facei;
                    }
                }
            }
            
            // 并行比较：找出全场最近距离
            if (Pstream::parRun())
            {
                scalar globalMinDist = minDist;
                reduce(globalMinDist, minOp<scalar>());
                
                // 只有距离足够接近全局最小值的处理器才保留 minFaceID
                if (mag(minDist - globalMinDist) > 1e-5) 
                {
                    minFaceID = -1; 
                }
            }

            // 判定触发条件
            if (minFaceID != -1 && TpatchPtr)
            {
                const fvPatchScalarField& Tpatch = *TpatchPtr;
                scalar faceT = Tpatch[minFaceID];
                
                bool tempReady = (faceT > site.Tactivate);
                bool timeReady = (currentTime - site.lastNucleationTime > site.waitTime);

                if (tempReady && timeReady)
                {
                    localTriggers[i] = true;
                    localBubbleCenters[i] = faceCenters[minFaceID];
                }
            }
        }
    }

    // =========================================================
    // 5. 全局状态同步 (Global Sync)
    // =========================================================
    if (Pstream::parRun())
    {
        // 使用 labelField 进行归约，这是最稳健的并行方法
        labelField triggerField(mySites_.size(), 0);
        forAll(localTriggers, i) 
        {
            if (localTriggers[i]) triggerField[i] = 1;
        }

        // 全局求和：只要有一个处理器触发，sum > 0
        reduce(triggerField, sumOp<labelField>());
        
        forAll(triggerField, i) 
        {
            if (triggerField[i] > 0) localTriggers[i] = true;
        }
        
        // 同步气泡真实的物理中心 (Sum Reduction)
        forAll(mySites_, i) 
        {
            if (localTriggers[i]) 
            {
                reduce(localBubbleCenters[i], sumOp<point>());
            }
        }
    }

    // =========================================================
    // 6. 执行生成与生长
    // =========================================================
    bool anyAction = false;

    forAll(mySites_, i)
    {
        NucleationSite& site = mySites_[i];

        // 6.1 激活新气泡
        if (localTriggers[i] && !site.isActive)
        {
            Info<< "BubbleGenerator: Nucleating at site " << i 
                << " (Center: " << localBubbleCenters[i] << ")" << endl;
            
            site.isActive = true;
            site.lifeTimer = 0;
            site.lastNucleationTime = currentTime;
            site.location = localBubbleCenters[i]; // 更新为实际面心坐标
        }

        // 6.2 气泡生长 (挖洞)
        if (site.isActive)
        {
            if (site.lifeTimer < maxLifeSteps_)
            {
                forAll(alpha, cI)
                {
                    scalar dist = mag(C[cI] - site.location);
                    if (dist < site.radius)
                    {
                        alpha[cI] = 0.0; // 设为气相
                        
                        // 重置温度为饱和温度 (消除过热液体导致的爆炸源项)
                        if (T[cI] > Tsat) 
                        {
                            T[cI] = Tsat;
                        }
                    }
                }
                site.lifeTimer++;
                anyAction = true;
            }
            else
            {
                site.isActive = false; // 寿命结束，停止保护
            }
        }
    }

    // 7. 修正边界与时间步
    if (anyAction)
    {
        alpha.correctBoundaryConditions();
        T.correctBoundaryConditions();
    }
    
    adjustDT(anyAction & fixDT_);

    return true;

    /* ---- Version 1.0 ----
       bool nucleated(false);
    const fvMesh& mesh = mesh_;
    const volVectorField& C = mesh.C();
    const scalar currentTime = runTime_.value();
    
    // 获取场指针
    volScalarField* ptrAlpha = const_cast<volScalarField*>(
        mesh.findObject<volScalarField>(alphaName_)
    );
    const volScalarField* ptrT = mesh.findObject<volScalarField>(TName_);

    if (!ptrAlpha || !ptrT)
    {
        // 如果找不到场，每步都会警告
        if (debug) 
        {
            Pout<< "BubbleGenerator Error: Fields " << alphaName_ 
                << " or " << TName_ << " not found!" << endl;
        }
        return false;
    }

    volScalarField& alpha = *ptrAlpha;
    const volScalarField& T = *ptrT;

    forAll(mySites_, i)
    {
        NucleationSite& site = mySites_[i];
        
        if (!site.isActive)
        {
            // 调试：打印搜索动作
            label cellID = mesh.findCell(site.location);

            // 只有找到该点的进程才会进入此逻辑
            if (cellID != -1)
            {
                scalar localT = T[cellID];
                scalar timePassed = currentTime - site.lastNucleationTime;
                bool tempReady = (localT > site.Tactivate);
                bool timeReady = (timePassed > site.waitTime);

                // --- 核心调试输出 ---
                Pout<< "Checking Site " << i << ": "
                    << "CellID=" << cellID 
                    << ", LocalT=" << localT 
                    << ", TargetT=" << site.Tactivate
                    << ", TimePassed=" << timePassed 
                    << ", WaitTime=" << site.waitTime << endl;

                if (tempReady && timeReady)
                {
                    Info<< ">>> NUCLEATING at site " << i 
                        << " (T=" << localT << " > " << site.Tactivate << ")" << endl;
                    site.isActive = true;
                    site.lifeTimer = 0;
                    site.lastNucleationTime = currentTime;
                    nucleated = true;
                }
                else if (tempReady && !timeReady)
                {
                    Pout<< "Site " << i << ": Temp OK, but waiting for cool-down." << endl;
                }
            }
        }

        // 执行挖洞逻辑
        if (site.isActive && site.lifeTimer < maxLifeSteps_)
        {
            Pout<< "Site " << i << ": Generating bubble, lifeTimer=" << site.lifeTimer << endl;
            forAll(alpha, cI)
            {
                if (mag(C[cI] - site.location) < site.radius)
                {
                    alpha[cI] = 0.0;
                }
            }
            site.lifeTimer++;
            nucleated = true;

            if (site.lifeTimer >= maxLifeSteps_)
            {
                site.isActive = false;
                Pout<< "Site " << i << ": Initialization finished." << endl;
            }
        }
    }

    bool globalNucleated = nucleated;
    if (Pstream::parRun())
    {
        // 使用 maxOp 代替 orOp，效果一样且兼容性更好
        // 逻辑：只要有一个进程是 true (1)，max 结果就是 true (1)
        reduce(globalNucleated, maxOp<bool>());
    }

    // 使用全局状态来决定是否调整时间步
    adjustDT(globalNucleated & fixDT_);

    return true;
    */

    /* --- code from Magnini et al. (2021) ignored ---
    bool nucleated(false);

    const vectorField& sites = location_->sites();

    forAll(sites,siteId)
    {
        if(activation_->isActive(sites[siteId]))
        {
            nucleated =
            (
                nucleated
                ||
                nucleation_->nucleate(sites[siteId])
            );
        }
    }

    adjustDT(nucleated & fixDT_);

    return true;
    */
}

bool Foam::functionObjects::bubbleGenerator::end()
{
    return true;
}


bool Foam::functionObjects::bubbleGenerator::write()
{
    writeState();
    return true;
}

void Foam::functionObjects::bubbleGenerator::adjustDT(bool adjust) const
{
    Time& runTime = const_cast<Time&>(runTime_);

    if(adjust)
    {
        Info<< "BubbleGenerator : Setting deltaT to " << newDT_
            << " s" << endl;
        runTime.setDeltaT(newDT_);
    }
}

void Foam::functionObjects::bubbleGenerator::writeState() const
{
    if (Pstream::master())
    {
        IOobject io
        (
            "nucleationState",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        );

        IOdictionary stateDict(io);

        forAll(mySites_, i)
        {
            const NucleationSite& site = mySites_[i];
            
            dictionary siteDict;
            siteDict.add("isActive", site.isActive);
            siteDict.add("lifeTimer", site.lifeTimer);
            siteDict.add("lastNucleationTime", site.lastNucleationTime);
            
            // [修正]
            // 1. 先拼接成 string，再显式构造成 word
            // 2. 然后再传给 add 函数
            word siteKey("site" + Foam::name(i));
            
            stateDict.add(siteKey, siteDict);
        }

        stateDict.regIOobject::write();
    }
}

void Foam::functionObjects::bubbleGenerator::readState()
{
    // 定义 IOobject，尝试读取
    IOobject io
    (
        "nucleationState",
        mesh_.time().timeName(),
        mesh_,
        IOobject::READ_IF_PRESENT, // 关键：如果文件存在就读，不存在（第一次计算）就跳过
        IOobject::NO_WRITE
    );

    // 检查文件头是否正常
    if (io.typeHeaderOk<IOdictionary>(true))
    {
        Info<< "BubbleGenerator: Restoring nucleation state from " << io.name() << " ..." << endl;
        
        IOdictionary stateDict(io);

        // 遍历所有站点
        forAll(mySites_, i)
        {
            NucleationSite& site = mySites_[i];
            word entryName = word("site") + Foam::name(i);

            // 检查字典里有没有这个站点的数据
            if (stateDict.found(entryName))
            {
                const dictionary& siteDict = stateDict.subDict(entryName);
                
                // 恢复状态
                siteDict.readEntry("isActive", site.isActive);
                siteDict.readEntry("lifeTimer", site.lifeTimer);
                siteDict.readEntry("lastNucleationTime", site.lastNucleationTime);
            }
        }
    }
    else
    {
        // 第一次计算，没有状态文件，保持默认初始化
        if (debug)
        {
            Info<< "BubbleGenerator: No previous state file found. Starting fresh." << endl;
        }
    }
}

// ************************************************************************* //
