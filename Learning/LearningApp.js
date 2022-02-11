covidsafe.viewSummaryPage = function (req, res, next) {
    const error = { location: "survey.covidsafe.viewSummaryPage" };

    const isStudent = true;

    const dbPerson = (isStudent) ? dbProfiles : dbFacstaff;
    const dataObj = Object.assign({}, req.session, { config });

    const date = serverDutyNight();
    const target = { date, isStudent };
    dbCovidsafe.get(target)
        .then(responses => {

            // create dictionary of responses by userID
            const responseDict = responses.reduce((acc, response) => {
                acc[response._id] = response;
                return acc;
            }, {});

            Object.assign(dataObj, target, { responseDict });
            const target = {};
            return dbPerson.get(target);
        })
        .then(population => {
            Object.assign(dataObj, { population });

            const responseDict = dataObj.responseDict;
            const missing = population.filter(person => {
                const response = responseDict[person.userID];
                return !response && person.userID > 0;
            })
            Object.assign(dataObj, { missing });

            const redList = population.filter(person => {
                const response = responseDict[person.userID];
                const isRed = response && !response.isCleared;
                if (isRed) {
                    person.raw = response.raw;
                }
                return isRed;
            })
            Object.assign(dataObj, { redList });

            const schema = new Question();
            const target = {};
            return dbSurveys.get("covidsafe", schema, target);
        })
        .then(questions => {

            // create dictionary of questions by id
            const questionDict = questions.reduce((acc, question) => {
                acc[question._id] = question.label;
                return acc;
            }, {});
            Object.assign(dataObj, { questionDict });

            res.status(201).render('pages/reportCovidsafe', dataObj);
        })
        .catch(err => {
            error.message = err.message || err.toString();
            return reportError(error, res, req.xhr);
        });
}