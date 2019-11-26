const yaml = require('js-yaml');
const fs = require('fs');

const args = process.argv;

keywords = args[2];
outfile = args[3];

try {
    const config = yaml.safeLoad(fs.readFileSync(keywords, 'utf8'));
    var results = {}
    for (var key in config) {
        results[key] = config[key].join("|")
    }
    fs.writeFileSync(outfile, JSON.stringify(results, null, 4));
} catch (e) {
    console.log(e);
}
