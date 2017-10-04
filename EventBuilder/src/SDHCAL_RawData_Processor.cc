#include "SDHCAL_RawData_Processor.h"
#include <iostream>
#include <iomanip>
#include <assert.h>

#include <EVENT/LCCollection.h>
#include <EVENT/LCGenericObject.h>
#include <UTIL/LCTOOLS.h>
#include "IMPL/LCCollectionVec.h"
#include "IMPL/LCFlagImpl.h"
#include "IMPL/RawCalorimeterHitImpl.h"

//#include <IMPL/LCFlagImpl.h>


//#include <IMPL/LCCollectionVec.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

#include "DIFUnpacker.h"


using namespace lcio ;
using namespace marlin ;



void SDHCAL_buffer::printBuffer(unsigned int start, unsigned int stop,std::ostream& flux)
{
  flux << std::hex;
  for (unsigned int k=start; k<stop; k++)
    flux << (unsigned int)(first[k]) << " - ";
  flux << std::dec <<  std::endl;
}




SDHCAL_RawBuffer_Navigator::SDHCAL_RawBuffer_Navigator(SDHCAL_buffer b) :
  _buffer(b),_SCbuffer(0,0)
{
  _DIFstartIndex=DIFUnpacker::getStartOfDIF(_buffer.buffer(),_buffer.getsize(),92);
  _theDIFPtr=NULL;
  _badSCdata=false;
}

DIFPtr* SDHCAL_RawBuffer_Navigator::getDIFPtr()
{
  if (NULL==_theDIFPtr) _theDIFPtr=new DIFPtr(getDIFBufferStart(),getDIFBufferSize());
  return _theDIFPtr;
}

uint32_t SDHCAL_RawBuffer_Navigator::getDIF_CRC()
{
  uint32_t i=getEndOfDIFData();
  uint32_t ret=0;
  ret |= ( (_buffer.buffer()[i-2])<<8 );
  ret |= _buffer.buffer()[i-1];
  return ret;
}

void SDHCAL_RawBuffer_Navigator::setSCBuffer()
{
  if (! hasSlowControlData() ) return;
  if (0 != _SCbuffer.getsize() ) return; //deja fait
  if (_badSCdata) return;
  _SCbuffer.first=&(getDIFBufferStart()[getEndOfDIFData()]);
  //compute Slow Control size
  uint32_t maxsize=_buffer.getsize()-_DIFstartIndex-getEndOfDIFData()+1; // should I +1 here ?
  uint32_t k=1; //SC Header
  uint32_t dif_ID=_SCbuffer.first[1];
  uint32_t chipSize=_SCbuffer.first[3];
    while ( (dif_ID != 0xa1 && _SCbuffer.first[k] != 0xa1 && k <maxsize) ||
	    (dif_ID == 0xa1 && _SCbuffer.first[k+2]==chipSize && k<maxsize) )
    {
      k+=2; //DIF ID + ASIC Header
      uint32_t scsize=_SCbuffer.first[k];
      if (scsize != 74 && scsize != 109)
	{
	  std::cout << "PROBLEM WITH SC SIZE " << scsize << std::endl;
	  k=0;
	  _badSCdata=true;
	  break;
	}
      k++; //skip size bit
      k+=scsize; // skip the data
    }
  if (_SCbuffer.first[k] == 0xa1 && !_badSCdata ) _SCbuffer.second=k+1; //add the trailer
  else
    {
      _badSCdata=true;
      std::cout << "PROBLEM SC TRAILER NOT FOUND " << std::endl;
    }
}


SDHCAL_buffer SDHCAL_RawBuffer_Navigator::getEndOfAllData()
{
  setSCBuffer();
  if (hasSlowControlData() && !_badSCdata)
    {
      return SDHCAL_buffer( &(_SCbuffer.buffer()[_SCbuffer.getsize()]), getSizeAfterDIFPtr()-3-_SCbuffer.getsize() );
    }
  else
    return SDHCAL_buffer( &(getDIFBufferStart()[getEndOfDIFData()]), getSizeAfterDIFPtr()-3 ); //remove the 2 bytes for CRC and the DIF trailer
}





SDHCAL_RawData_Processor aSDHCAL_RawData_Processor ;


SDHCAL_RawData_Processor::SDHCAL_RawData_Processor() : Processor("SDHCAL_RawData_Processor") {

    // modify processor description
    _description = "SDHCAL_RawData_Processor converts SDHCAL Raw Data into collection of RawCalorimeterHit" ;


    // register steering parameters: name, description, class-variable, default value

    bool debugDefault=false;
    _debugMode=false;
    registerProcessorParameter( "DebugMode",
				"Turn ON/OFF debug mode : Warning Debug mode uses assert and may crash the application",
				_debugMode,
				debugDefault);

   registerInputCollection( LCIO::LCGENERICOBJECT,
			      "XDAQCollectionName" ,
			      "XDAQ produced collection name" ,
			      _XDAQCollectionNames,
			     std::string("RU_XDAQ"));



    registerOutputCollection( LCIO::RAWCALORIMETERHIT,
			      "OutputRawCaloHitCollectionName" ,
			      "Name of output collection containing raw calorimeter hits"  ,
			      _RawHitCollectionName ,
			      std::string("DHCALRawHits")
			      );


    _nrun=_nevt=_nWrongObj=_nProcessedObject=_hasSlowControl=_hasBadSlowControl=0;
}



void SDHCAL_RawData_Processor::init() {

  //streamlog_out(DEBUG) << "   init called  " << std::endl ;

  // usually a good idea to
  printParameters() ;


}


void SDHCAL_RawData_Processor::processRunHeader( LCRunHeader* run) {
  _nrun++;
  LCTOOLS::dumpRunHeader(run);
}


void SDHCAL_RawData_Processor::processEvent( LCEvent * evt ) {

  _nevt++;

  IMPL::LCCollectionVec *RawVec=new IMPL::LCCollectionVec(LCIO::RAWCALORIMETERHIT) ;

  //Prepare a flag to tag data type in RawVec (dit les types de data qu'on va enregistrer?)
  IMPL::LCFlagImpl chFlag(0) ;
  EVENT::LCIO bitinfo;
  chFlag.setBit(bitinfo.RCHBIT_LONG ) ;                    // raw calorimeter data -> format long //(sert a qq chose?)
  chFlag.setBit(bitinfo.RCHBIT_BARREL ) ;                  // barrel
  chFlag.setBit(bitinfo.RCHBIT_ID1 ) ;                     // cell ID
  chFlag.setBit(bitinfo.RCHBIT_TIME ) ;                    // timestamp
  RawVec->setFlag( chFlag.getFlag()  ) ;



  try{
    LCCollection* col = evt->getCollection( _XDAQCollectionNames );
    int nElement=col->getNumberOfElements();
    _CollectionSizeCounter[nElement]++;
    for (int iel=0; iel<nElement; iel++)
      {
	LCGenericObject* obj=dynamic_cast<LCGenericObject*>(col->getElementAt(iel));
	if (NULL==obj)
	  {
	    _nWrongObj++;
	    continue;
	  }
	_nProcessedObject++;
	LMGeneric *lmobj=(LMGeneric *) obj;
	unsigned char* debug_variable_1=lmobj->getSDHCALBuffer().endOfBuffer();
	SDHCAL_RawBuffer_Navigator bufferNavigator(lmobj->getSDHCALBuffer());
	unsigned char* debug_variable_2=bufferNavigator.getDIFBuffer().endOfBuffer();
	streamlog_out(DEBUG) << "DIF BUFFER END " << (unsigned int *) debug_variable_1 << " " << (unsigned int *) debug_variable_2 << std::endl;
	if (_debugMode) assert (debug_variable_1 == debug_variable_2);

	uint32_t idstart=bufferNavigator.getStartOfDIF();
	streamlog_out(DEBUG) << "DIF header starts at " << idstart
			     << " and is equal to " << std::hex << (unsigned int) lmobj->getCharBuffer()[idstart] << std::dec << std::endl;
	if(idstart==0 && streamlog::out.write< streamlog::DEBUG >() ) lmobj->getSDHCALBuffer().printBuffer( 0 , streamlog::out() );
	_DIFStarter[idstart]++;
	if (!  bufferNavigator.validBuffer() ) continue;
	if (_debugMode)
	  {
	    //for PS and SPS data in april-may 2012
	    if (idstart!=0)
	      {
		assert(idstart==94);  // DIF global header is at byte 94
	      }
	    else
	      assert(lmobj->nBytes()==224); //Don't know why but for the moment will assume that
	  }

	DIFPtr *d=bufferNavigator.getDIFPtr();
	if (_debugMode) assert(d!=NULL);
	if (d!=NULL)
	  {
	    _DIFPtrValueAtReturnedPos[bufferNavigator.getDIFBufferStart()[d->getGetFramePtrReturn()] ]++;
	    if (_debugMode) assert( bufferNavigator.getDIFBufferStart()[d->getGetFramePtrReturn()]==0xa0 );
	  }
	_SizeAfterDIFPtr[bufferNavigator.getSizeAfterDIFPtr()]++;

	//createRawCalorimeterHit
	for (uint32_t i=0;i<d->getNumberOfFrames();i++)
	  {
	    for (uint32_t j=0;j<64;j++)
	      {
		if (!(d->getFrameLevel(i,j,0) || d->getFrameLevel(i,j,1))) continue; // skip empty pads
		//  std::cout <<" New hit "<<std::endl;
		unsigned long int ID0;
		ID0=(unsigned long int)(((unsigned short)d->getID())&0xFF);			//8 firsts bits: DIF Id
		ID0+=(unsigned long int)(((unsigned short)d->getFrameAsicHeader(i)<<8)&0xFF00);	//8 next bits:   Asic Id
		bitset<6> Channel(j);
		ID0+=(unsigned long int)((Channel.to_ulong()<<16)&0x3F0000);			//6 next bits:   Asic's Channel
		unsigned long BarrelEndcapModule=0;  //(40 barrel + 24 endcap) modules to be coded here  0 for testbeam (over 6 bits)
		ID0+=(unsigned long int)((BarrelEndcapModule<<22)&0xFC00000);
		unsigned long int ID1 = (unsigned long int)(d->getFrameBCID(i));
		bitset<2> ThStatus;
		ThStatus[0]=d->getFrameLevel(i,j,0);
		ThStatus[1]=d->getFrameLevel(i,j,1);
		//ThStatus[2]=isSynchronised; //I'm not computing the synchronisation here.

		IMPL::RawCalorimeterHitImpl *hit=new IMPL::RawCalorimeterHitImpl() ;
		hit->setCellID0((unsigned long int)ID0);
		hit->setCellID1(ID1);
		hit->setAmplitude(ThStatus.to_ulong());
		unsigned long int TTT = (unsigned long int)(d->getFrameTimeToTrigger(i));
		hit->setTimeStamp(TTT);			      		//Time stamp of this event from Run Begining
		RawVec->addElement(hit);
	      }//for (uint32_t j=0;j<64;j++)

	  }//for (uint32_t i=0;i<d->getNumberOfFrames();i++)

	//register Triggers'time : lots of limits.here ?

	lcio::IntVec trig(8);
	trig[0] = d->getDTC();
	trig[1] = d->getGTC();
	trig[2] = d->getBCID();
	trig[3] = d->getAbsoluteBCID()&0xFFFFFF;
	trig[4] = (d->getAbsoluteBCID()/(0xFFFFFF+1))&0xFFFFFF;
	trig[5] = d->getTASU1(); //what happen if no temperature ?
	trig[6] = d->getTASU2();
	trig[7] = d->getTDIF();

	std::stringstream ss("");
	ss<<"DIF"<<d->getID()<<"_Triggers";
	RawVec->parameters().setValues(ss.str(),trig);


	if (bufferNavigator.hasSlowControlData()) _hasSlowControl++;
	if (bufferNavigator.badSCData())  _hasBadSlowControl++;

	SDHCAL_buffer eod=bufferNavigator.getEndOfAllData();
	_SizeAfterAllData[eod.getsize()]++;
	unsigned char* debug_variable_3=eod.endOfBuffer();
	streamlog_out(DEBUG) << "END DATA BUFFER END " << (unsigned int *) debug_variable_1 << " " << (unsigned int *) debug_variable_3 << std::endl;
	if (_debugMode) assert (debug_variable_1 == debug_variable_3);
	streamlog_out(DEBUG) << "End of Data remaining stuff : ";
	if(streamlog::out.write< streamlog::DEBUG >() ) eod.printBuffer( 0 , streamlog::out() );

	int nonzeroCount=0;
	for (unsigned char* it=eod.buffer(); it != eod.endOfBuffer(); it++)
	  if (int(*it) !=0) nonzeroCount++;
	_NonZeroValusAtEndOfData[nonzeroCount]++;


      } //for (int iel=0; iel<nElement; iel++)

  }
  catch(DataNotAvailableException &e){
    streamlog_out(WARNING) << _XDAQCollectionNames << " collection not available" << std::endl;
  }


  evt->addCollection(RawVec,_RawHitCollectionName);
}



void SDHCAL_RawData_Processor::check( LCEvent * evt ) {
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void SDHCAL_RawData_Processor::printCounter(std::string description, std::map<int,int> &m)
{
  streamlog_out(MESSAGE) << " statistics for " << description << " : ";
  for (std::map<int,int>::iterator it=m.begin(); it != m.end(); ++it)
    {
      if (it != m.begin() ) streamlog_out(MESSAGE) << ",";
      streamlog_out(MESSAGE) << " [" << it->first << "]=" << it->second;
    }
  streamlog_out(MESSAGE) << std::endl;
}


void SDHCAL_RawData_Processor::end(){

  streamlog_out(MESSAGE) << "FINAL STATISTICS : " << std::endl;
  streamlog_out(MESSAGE) << _nrun << " runs for " << _nevt << " events." << std::endl;
  printCounter("Size of GenericObject collections",_CollectionSizeCounter);
  streamlog_out(MESSAGE) << "Number of badly extracted LCGenericObject : " <<  _nWrongObj << std::endl;
  streamlog_out(MESSAGE) << "Number of processed LCGenericObject : " <<  _nProcessedObject << std::endl;
  printCounter("Start of DIF header",_DIFStarter);
  printCounter("Value after DIF data are processed",_DIFPtrValueAtReturnedPos);
  printCounter("Size remaining in buffer after end of DIF data",  _SizeAfterDIFPtr );
  streamlog_out(MESSAGE) << "Number of Slow Control found " << _hasSlowControl << " out of which " << _hasBadSlowControl << " are bad" << std::endl;
  printCounter("Size remaining after all of data have been processed",_SizeAfterAllData );
  printCounter("Number on non zero values in end of data buffer",_NonZeroValusAtEndOfData);

}

