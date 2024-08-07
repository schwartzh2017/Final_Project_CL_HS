const fs = require('fs');
const { parse } = require('csv-parse');
const { stringify } = require('csv-stringify');
const { calculateRisk } = require('./riskCalculator');

const inputFile = 'input.csv';
const outputFile = 'output.csv';

const inputStream = fs.createReadStream(inputFile);
const outputStream = fs.createWriteStream(outputFile);

const parser = parse({ columns: true });
const stringifier = stringify({ header: true });

inputStream.pipe(parser);
stringifier.pipe(outputStream);

parser.on('data', (row) => {
  const age = parseInt(row.age);
  const race = parseInt(row.race);
  const relativeBC = parseInt(row.relativeBC);
  const biopsy = parseInt(row.biopsy);
  const density = parseInt(row.density);

  const result = calculateRisk(age, race, relativeBC, biopsy, density);

  stringifier.write({
    ...row,
    risk_5Y: result.risk_5Y,
    risk_10Y: result.risk_10Y,
    risk_5Y_average: result.risk_5Y_average,
    risk_10Y_average: result.risk_10Y_average
  });
});

parser.on('end', () => {
  stringifier.end();
  console.log('Risk calculation completed. Results written to output.csv');
});
