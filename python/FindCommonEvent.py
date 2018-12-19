#!/usr/bin/env python2

from __future__ import print_function
from ROOT import TCanvas, TMath, TH2F, TFile, TTree, TList, THStack, TObject, gDirectory, gPad, kRed
import json
import itertools
import platform

runNumber = 744221
maxAcquisition = 1000  # Cut on number of trig
timeDif = 500  # in clocks

# sdhcalCuts
maxCogDeviation = 50  # in mm, dif between cog of ecal/sdhcal hits
sdhcalCutX = [225, 401]
sdhcalCutY = [377, 553]
minHits = 20

# ecalCuts
minSlabs = 3

runOnLxplus = False  # Automatically adapt ilcsoft paths when on lxplus
if platform.node().find("lxplus") != -1:
    runOnLxplus = True

sdhcalDataPath = "/eos/user/a/apingaul/CALICE/Data/SPS_09_2018/Ecal/"
if runOnLxplus is True:
    ecalDataPath = "/eos/project/s/siw-ecal/TB2018-09/Common/ECAL/offset_twiki/Muon_200GeV/"
else:
    ecalDataPath = sdhcalDataPath

hcalFile = TFile.Open(sdhcalDataPath + 'TDHCAL_Ecal_{}.root'.format(runNumber), 'read')
hcalTreeName = 'sdhcal'
hcalSelectionCut = 'HitCogX>{} && HitCogX<{} && HitCogY>{} && HitCogY<{} && NHits>{} && TrigNum<{}'.format(
    sdhcalCutX[0], sdhcalCutX[1], sdhcalCutY[0], sdhcalCutY[1], minHits, maxAcquisition
    if maxAcquisition != 0 else 10000)

ecalFile = TFile.Open(ecalDataPath + '{}__build.root'.format(runNumber), 'read')
# ecalFile = TFile.Open(ecalDataPath + 'ECAL_200GeV_muons_build.root', 'read')
ecalTreeName = 'ecal'
ecalSelectionCut = 'nhit_slab>{} && spill<{}'.format(minSlabs, maxAcquisition if maxAcquisition != 0 else 10000)
hcalOutFile = TFile.Open(sdhcalDataPath + 'TDHCAL_Common_{}.root'.format(runNumber), 'recreate')
# sdhcal->Scan("TrigNum:EvtBcid:EvtRevBcid:HitCogX:HitCogY:NHits","HitCogX>225 && HitCogX<401 && HitCogY>377 && HitCogY<553","box", 2541, 0);
# ecal->Scan("spill:prev_bcid:bcid:next_bcid:nhit_slab", "nhit_slab>3", "", 12962, 0);
hTree = hcalFile.Get(hcalTreeName).CopyTree(hcalSelectionCut)
eTree = ecalFile.Get(ecalTreeName).CopyTree(ecalSelectionCut)

commonEvtList = []
commonRevEvtList = []

for evt in eTree:
    hSmallTree = hTree.CopyTree("TrigNum == {}".format(evt.spill))
    for hcalEvt in hSmallTree:
        assert (evt.spill == hcalEvt.TrigNum)
        # print ("spill: '{}' evt '{}' - trig: '{}' evt'{}'".format(evt.spill, evt.event, hcalEvt.TrigNum, hcalEvt.EvtNum))
        if abs(hcalEvt.EvtBcid - evt.bcid) < timeDif:
            commonEvtList.append({
                "ecal": {
                    "evtNum": evt.event,
                    "bcid": evt.bcid
                },
                "hcal": {
                    "evtNum": hcalEvt.EvtNum,
                    "bcid": hcalEvt.EvtBcid
                }
            })
            # print(evt.spill, hcalEvt.TrigNum, evt.event, hcalEvt.EvtNum, evt.prev_bcid, evt.bcid, evt.next_bcid, evt.nhit_slab)
            # print(evt.spill, hcalEvt.TrigNum, evt.event, hcalEvt.EvtNum, abs(evt.bcid - hcalEvt.EvtBcid))
        elif abs(hcalEvt.EvtRevBcid - evt.bcid) < timeDif:
            commonRevEvtList.append({
                "ecal": {
                    "evtNum": evt.event,
                    "bcid": evt.bcid
                },
                "hcal": {
                    "evtNum": hcalEvt.EvtNum,
                    "bcid": hcalEvt.EvtRevBcid
                }
            })
            # print(evt.spill, hcalEvt.TrigNum, evt.event, hcalEvt.EvtNum, abs(evt.bcid - hcalEvt.EvtRevBcid))

# print("common: ")
# print(json.dumps(commonEvtList, sort_keys=True, indent=2))

# print("\n RevCommon: ")
# print(json.dumps(commonRevEvtList, sort_keys=True, indent=2))

commonTree = TTree('commonTree', 'tree with potential common evts')
# commonRevTree = TTree('commonEvtRevTime', 'tree with potential common evts, with hcal reversed time')
eList = TList()
hList = TList()
for evt in commonRevEvtList:
    eEvtTree = eTree.CopyTree("event == {}".format(evt['ecal']['evtNum']))
    eEvtTree.Scan(
        "spill:event:bcid"
    )  # For some root magic reason if I dont scan the event variable here it is filled with some 'random' value, other variables are fine.......
    assert (evt['ecal']['evtNum'] == eEvtTree.event)
    assert (evt['ecal']['bcid'] == eEvtTree.bcid)
    # print(eEvtTree.spill, eEvtTree.event, eEvtTree.bcid)
    eList.Add(eEvtTree)

    hEvtTree = hTree.CopyTree("EvtNum == {}".format(evt['hcal']['evtNum']))
    hEvtTree.Scan(
        "TrigNum:EvtNum:EvtRevBcid:NHits"
    )  # For some root magic reason if I dont scan the EvtNum variable here it is filled with some 'random' value, other variables are fine.......
    assert (evt['hcal']['evtNum'] == hEvtTree.EvtNum)
    assert (evt['hcal']['bcid'] == hEvtTree.EvtRevBcid)
    # print(hEvtTree.TrigNum, hEvtTree.EvtNum, hEvtTree.EvtBcid)
    hList.Add(hEvtTree)

assert eList.LastIndex() != -1 and hList.LastIndex() != -1, "No common event Found"
eTreeOut = TTree.MergeTrees(eList)
eTreeOut.SetName('ecal')
eTreeOut.Write()
hTreeOut = TTree.MergeTrees(hList)
hTreeOut.SetName('sdhcal')
hTreeOut.Write()
commonTree.AddFriend(eTreeOut)
commonTree.AddFriend(hTreeOut)
commonTree.Write()

# for comEvt in eTreeOut:
# for eEvt, hEvt in itertools.izip(ecal, sdhcal):
i = 0
can = TCanvas("name{}".format(i), "title{}".format(i), 1800, 1000)
can.Divide(0, 2)
for eEvt, hEvt in itertools.izip(eTreeOut, hTreeOut):
    entry = eEvt.GetReadEntry()
    # commonTree.Scan("*", "Entry$=={}".format(entry))
    eEvtNum = commonTree.event
    hEvtNum = commonTree.EvtNum
    # print(commonTree.spill, commonTree.event, commonTree.EvtNum, commonTree.bcid, commonTree.EvtRevBcid,
    #       abs(commonTree.bcid - commonTree.EvtRevBcid))  # , commonTree.NHits)
    hCogX = commonTree.HitCogX
    hCogY = commonTree.HitCogY
    eCogX = 225 + 2.25 + 5.5 * TMath.Mean(len(commonTree.hit_x), commonTree.hit_x)
    eCogY = 377. + 2.25 + 5.5 * TMath.Mean(len(commonTree.hit_y), commonTree.hit_y)
    # eCogX = 313 + TMath.Mean(len(commonTree.hit_x), commonTree.hit_x)
    # eCogY = 465 + TMath.Mean(len(commonTree.hit_y), commonTree.hit_y)
    # print(hCogX, eCogX, abs(hCogX - eCogX), hCogY, eCogY, abs(hCogY - eCogY))

    if abs(eCogY - hCogY) < maxCogDeviation and abs(eCogX - hCogX) < maxCogDeviation:
        print(commonTree.spill, commonTree.event, commonTree.EvtNum, commonTree.bcid, commonTree.EvtRevBcid,
              abs(commonTree.bcid - commonTree.EvtRevBcid))  # , commonTree.NHits)

        deltaCogX = int(abs(hCogX - eCogX))
        deltaCogY = int(abs(hCogY - eCogY))
        deltaBcid = abs(commonTree.bcid - commonTree.EvtRevBcid)
        print("DeltaCogx: ", deltaCogX, " - DeltaCogy: ", deltaCogY)

        can.cd(1)
        commonTree.Draw("225+5.5*hit_x+2.25:hit_z*15>>xecal", "Entry$ == {}".format(entry), "goff")
        # commonTree.Draw("313+hit_x:hit_z*15>>xecal", "Entry$ == {}".format(entry), "goff")
        commonTree.Draw("HitX:HitZ>>xhcal", "Entry$ == {}".format(entry), "goff")
        hs = THStack("hs{}".format(i),
                     "TopView - #Delta CogX = {}mm - #Delta CogY = {}mm - #Delta Bcid = {}clocks;z (mm);x (mm)".format(
                         deltaCogX, deltaCogY, deltaBcid))
        xecal = gDirectory.Get("xecal")
        xhcal = gDirectory.Get("xhcal")
        xecal.SetLineColor(kRed)
        xecal.GetYaxis().SetLimits(200, 450)
        xecal.GetXaxis().SetLimits(0, 1250)
        xhcal.GetYaxis().SetLimits(200, 450)
        xhcal.GetXaxis().SetLimits(0, 1250)
        can.Update()
        hs.Add(xecal)
        hs.Add(xhcal)
        hs.Draw("box")
        can.Update()

        can.cd(2)
        commonTree.Draw("377+5.5*hit_y+2.25:hit_z*15>>yecal", "Entry$ == {}".format(entry), "goff")
        # commonTree.Draw("465+hit_y:hit_z*15>>yecal", "Entry$ == {}".format(entry), "goff")
        commonTree.Draw("HitY:HitZ>>yhcal", "Entry$ == {}".format(entry), "goff")
        hs1 = THStack("hs1{}".format(i), "SideView;z (mm);y (mm)")
        yecal = gDirectory.Get("yecal")
        yhcal = gDirectory.Get("yhcal")
        yecal.SetLineColor(kRed)
        yecal.GetYaxis().SetLimits(350, 600)
        yecal.GetXaxis().SetLimits(0, 1250)
        yhcal.GetYaxis().SetLimits(350, 600)
        yhcal.GetXaxis().SetLimits(0, 1250)
        can.Update()
        hs1.Add(yecal)
        hs1.Add(yhcal)
        hs1.Draw("box")
        can.Update()

        i += 1

        # meanX = TMath.Mean(commonTree.hit_x.size())
        raw_input("Press Enter to continue...")

hcalOutFile.Close()
