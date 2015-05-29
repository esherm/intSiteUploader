-- MySQL dump 10.13  Distrib 5.5.43, for debian-linux-gnu (x86_64)
--
-- Host: microb98.med.upenn.edu    Database: intsitesdev
-- ------------------------------------------------------
-- Server version	5.5.5-10.0.15-MariaDB

/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8 */;
/*!40103 SET @OLD_TIME_ZONE=@@TIME_ZONE */;
/*!40103 SET TIME_ZONE='+00:00' */;
/*!40014 SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0 */;
/*!40014 SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0 */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='NO_AUTO_VALUE_ON_ZERO' */;
/*!40111 SET @OLD_SQL_NOTES=@@SQL_NOTES, SQL_NOTES=0 */;

--
-- Table structure for table `multihitlengths`
--

DROP TABLE IF EXISTS `multihitlengths`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `multihitlengths` (
  `multihitClusterID` int(11) unsigned NOT NULL,
  `length` int(11) unsigned NOT NULL,
  `count` int(8) unsigned NOT NULL DEFAULT '1',
  KEY `multihitClusterID` (`multihitClusterID`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `pcrbreakpoints`
--

DROP TABLE IF EXISTS `pcrbreakpoints`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `pcrbreakpoints` (
  `siteID` int(11) unsigned NOT NULL,
  `breakpoint` int(11) unsigned NOT NULL,
  `count` int(8) unsigned NOT NULL DEFAULT '1',
  PRIMARY KEY (`siteID`,`breakpoint`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `samples`
--

DROP TABLE IF EXISTS `samples`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `samples` (
  `sampleID` int(11) unsigned NOT NULL,
  `sampleName` varchar(255) NOT NULL DEFAULT '',
  `refGenome` varchar(10) NOT NULL DEFAULT '',
  `gender` char(1) NOT NULL DEFAULT '',
  PRIMARY KEY (`sampleID`),
  UNIQUE KEY `sampleName` (`sampleName`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `sites`
--

DROP TABLE IF EXISTS `sites`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `sites` (
  `siteID` int(11) unsigned NOT NULL,
  `sampleID` int(11) unsigned NOT NULL,
  `position` int(11) unsigned NOT NULL,
  `chr` varchar(20) NOT NULL DEFAULT '',
  `strand` char(1) NOT NULL DEFAULT '',
  `multihitClusterID` int(11) DEFAULT NULL,
  PRIMARY KEY (`siteID`),
  KEY `sampleID` (`sampleID`),
  KEY `multihitClusterID` (`multihitClusterID`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping routines for database 'intsitesdev'
--
/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;

/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40014 SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;

-- Dump completed on 2015-05-29 12:00:03
