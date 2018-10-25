#!/usr/bin/env python2

from __future__ import print_function
from ROOT import TCanvas, TH2F, TFile, TTree, TList, THStack, TObject, gDirectory, gPad
import json
import itertools


runNumber = 744211
maxAcquisition = 0
timeDif = 500  # in clocks

hcalFile = TFile.Open('/eos/user/a/apingaul/CALICE/Data/SPS_09_2018/Trivent/TDHCAL_{}.root'.format(runNumber), 'read')
hcalTreeName = 'sdhcal'
hcalSelectionCut = 'HitCogX>225 && HitCogX<401 && HitCogY>377 && HitCogY<553 && NHits>20 && TrigNum<{}'.format(maxAcquisition if maxAcquisition != 0 else 10000)

ecalFile = TFile.Open('/eos/user/a/apingaul/CALICE/Data/SPS_09_2018/Ecal/{}__build.root'.format(runNumber), 'read')
ecalTreeName = 'ecal'
ecalSelectionCut = 'nhit_slab>3 && spill<{}'.format(maxAcquisition if maxAcquisition != 0 else 10000)


hcalOutFile = TFile.Open('/eos/user/a/apingaul/CALICE/Data/SPS_09_2018/Trivent/TDHCAL_Common_{}.root'.format(runNumber), 'recreate')
# sdhcal->Scan("TrigNum:EvtBcid:EvtRevBcid:HitCogX:HitCogY:NHits","HitCogX>225 && HitCogX<401 && HitCogY>377 && HitCogY<553","box", 2541, 0);
# ecal->Scan("spill:prev_bcid:bcid:next_bcid:nhit_slab", "nhit_slab>3", "", 12962, 0);
hTree = hcalFile.Get(hcalTreeName).CopyTree(hcalSelectionCut)
eTree = ecalFile.Get(ecalTreeName).CopyTree(ecalSelectionCut)

commonEvtList = []
commonRevEvtList = []

for evt in eTree:
    hSmallTree = hTree.CopyTree("TrigNum == {}".format(evt.spill))
    for hcalEvt in hSmallTree:
        assert(evt.spill == hcalEvt.TrigNum)
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

print("common: ")
print(json.dumps(commonEvtList, sort_keys=True, indent=2))

print("\n RevCommon: ")
print(json.dumps(commonRevEvtList, sort_keys=True, indent=2))

commonTree = TTree('commonTree', 'tree with potential common evts')
# commonRevTree = TTree('commonEvtRevTime', 'tree with potential common evts, with hcal reversed time')
eList = TList()
hList = TList()
for evt in commonRevEvtList:
    eEvtTree = eTree.CopyTree("event == {}".format(evt['ecal']['evtNum']))
    eEvtTree.Scan("spill:event:bcid")  # For some root magic reason if I dont scan the event variable here it is filled with some 'random' value, other variables are fine.......
    assert (evt['ecal']['evtNum'] == eEvtTree.event)
    assert (evt['ecal']['bcid'] == eEvtTree.bcid)
    # print(eEvtTree.spill, eEvtTree.event, eEvtTree.bcid)
    eList.Add(eEvtTree)

    hEvtTree = hTree.CopyTree("EvtNum == {}".format(evt['hcal']['evtNum']))
    hEvtTree.Scan("TrigNum:EvtNum:EvtRevBcid")  # For some root magic reason if I dont scan the EvtNum variable here it is filled with some 'random' value, other variables are fine.......
    assert (evt['hcal']['evtNum'] == hEvtTree.EvtNum)
    assert (evt['hcal']['bcid'] == hEvtTree.EvtRevBcid)
    # print(hEvtTree.TrigNum, hEvtTree.EvtNum, hEvtTree.EvtBcid)
    hList.Add(hEvtTree)

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
    print(commonTree.spill, commonTree.event, commonTree.EvtNum, commonTree.bcid, commonTree.EvtRevBcid, abs(commonTree.bcid - commonTree.EvtRevBcid))

    can.cd(1)
    commonTree.Draw("313+hit_x:hit_z*15>>xecal", "Entry$ == {}".format(entry), "goff")
    commonTree.Draw("HitX:HitZ>>xhcal", "Entry$ == {}".format(entry), "goff")
    hs = THStack("hs{}".format(i), "hsta{}".format(i))
    xecal = gDirectory.Get("xecal")
    xhcal = gDirectory.Get("xhcal")
    xecal.GetYaxis().SetLimits(0, 600)
    xecal.GetXaxis().SetLimits(0, 1200)
    can.Update()
    hs.Add(xecal)
    hs.Add(xhcal)
    hs.Draw("box nostack")
    can.Update()

    can.cd(2)
    commonTree.Draw("465+hit_y:hit_z*15>>yecal", "Entry$ == {}".format(entry), "goff")
    commonTree.Draw("HitY:HitZ>>yhcal", "Entry$ == {}".format(entry), "goff")
    hs1 = THStack("hs1{}".format(i), "hsta1{}".format(i))
    yecal = gDirectory.Get("yecal")
    yhcal = gDirectory.Get("yhcal")
    yecal.GetYaxis().SetLimits(0, 600)
    yecal.GetXaxis().SetLimits(0, 1200)
    can.Update()
    hs1.Add(yecal)
    hs1.Add(yhcal)
    hs1.Draw("box nostack")
    can.Update()

    i += 1
    raw_input("Press Enter to continue...")


hcalOutFile.Close()
