import open from "open";
import fs from "fs";
import { parse } from "csv-parse";

const url = (process.argv[2] || "../controllare.csv").replace("\\", "/");

fs.createReadStream(url)
    .pipe(parse({ delimiter: ";", from_line: 2 }))
    .on("data", async function (row) {
        // console.log("\n" + row[0]);
        if (row[1] != "")
            await open((row[1][0] == "h" ? "" : "https://doi.org/") + row[1]);
    }).on("end", function () {
        console.log("finished");
    }).on("error", function (error) {
        console.log(error.message);
    });