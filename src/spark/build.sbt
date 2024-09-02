name := "iqlib-spark"
version := "1.0"
scalaVersion := "2.13.0"
offline := true

scalacOptions += "-target:jvm-1.7"
javacOptions ++= Seq("-source", "1.7")



sparkComponents := Seq("core", "sql","graphx")
sparkVersion := "3.5.0"

// libraryDependencies ++= Seq(
//   "org.scala-lang.modules" % "scala-xml_2.11" % "1.0.5",
//   "org.apache.spark" %% "spark-core" % sparkVersion.value,
//   "org.apache.spark" %% "spark-sql" % sparkVersion.value,
//   "org.apache.spark" %% "spark-graphx" % sparkVersion.value
// )

libraryDependencies ++= Seq(
  // "org.scala-lang.modules" %% "scala-parser-combinators" % "1.1.0",
  // "org.scala-lang.modules" % "scala-xml_2.11" % "1.0.5",
  "org.apache.spark" %% "spark-core" % sparkVersion.value,
  "org.apache.spark" %% "spark-sql" % sparkVersion.value,
  "org.apache.spark" %% "spark-graphx" % sparkVersion.value
  // "org.locationtech.geotrellis" %% "geotrellis-raster" % "1.0.0",
  // "org.locationtech.geotrellis" %% "geotrellis-spark" % "1.0.0",
  // "org.vegas-viz" %% "vegas-spark" % "0.3.11"
)






spIgnoreProvided := true
