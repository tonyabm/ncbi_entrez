from flask import Flask, request, jsonify, render_template
from Bio import Entrez

app = Flask(__name__)

# 设置Entrez的电子邮件（NCBI要求）
Entrez.email = "tony.liu.abm@gmail.com"


# 首页路由，显示表单
@app.route('/')
def index():
    return render_template('index.html')


@app.route('/search_accession', methods=['POST'])
def search_accession():
    try:
        accession = request.form.get('accession')
        gene_name = request.form.get('gene_name')
        gene_id   = request.form.get('gene_id')

        if gene_id != '':
            handle = Entrez.efetch(db="gene", id=gene_id, rettype="gb", retmode="text")
            data = handle.read()
        elif gene_name != '':
            handle = Entrez.efetch(db="gene", id=gene_name, rettype="gb", retmode="text")
            data = handle.read()
        elif accession != '':
            handle = Entrez.efetch(db="nuccore", id=accession, rettype="gb", retmode="text")
            data = handle.read()
        handle.close()
        return jsonify({'status': 'ok', 'data': data})
    except:
        return jsonify({'status': 'error', 'data': 'Failed to retrieve data'})
    
@app.route('/pubmed_search', methods=['POST'])
def pubmed_search():

    # import dependencies
    import requests

    s = request.form.get('pubmed')
    # this is the base url of NCBI's Entrez API
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

    # accessing NCBI database
    database = "pubmed"

    # build query url that is used to retrieve search results
    query_url = (
            base_url + "esearch.fcgi?db=" + database + "&term=" + s + "&retmode=json"
    )

    # request in response for the search information in JSON format
    response = requests.get(query_url).json()

    # idList itself is a list, create a string with the same name
    idList = response["esearchresult"]["idlist"]

    allArticles = []

    # set condition if no articles are found
    if len(idList) == 0:
        print("No articles found! Please try another search term.")

    # condition when article(s) are found
    else:
        # for loop to iterate the first five articles in idList
        for article in idList[0:5]:

            # count articles
            total = len(idList[0:5])

            # find the index of each article
            i = idList.index(article)

            # message to show the article search
            if i == 0:
                print("------------------------")
                print("Beginning article search")
                print("------------------------")

            # create article count. python index starts at 0 so we need to add 1
            article_count = int(i) + 1

            # count number of each article
            print(f"Retrieving articles {article_count} of {total}")

            # base url to retrieve article summary
            url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?"

            # generate endpoint url from our search to retrieve information for each article
            search_url = url + "db=" + database + "&id=" + article + "&retmode=json"

            # request article information, formatted as json
            search_response = requests.get(search_url).json()

            # create empty list for storing multiple authors per article
            author_list = []

            # try and except loop to keep the loop running
            try:

                # data parsing from our search response
                pubmed_id = search_response["result"][article]["uid"]
                title = search_response["result"][article]["title"]
                authors = search_response["result"][article]["authors"]
                journal = search_response["result"][article]["source"]
                pub_date = search_response["result"][article]["pubdate"]
                volume = search_response["result"][article]["volume"]
                issue = search_response["result"][article]["issue"]
                pages = search_response["result"][article]["pages"]
                doi = search_response["result"][article]["elocationid"]

                # append values to list
                for i in authors:
                    all_authors = i["name"]
                    author_list.append(all_authors)

                # replace [,], and ' symbols from each element in author list
                names = (
                    str(author_list).replace("'", "").replace("[", "").replace("]", "")
                )

                # replace italization in title with <i></i>
                corrected_title = title.replace("&lt;i&gt;", "<i>").replace(
                    "&lt;/i&gt;", "</i>"
                )

                allArticles.append(
                    {
                        "pubmed_id": pubmed_id,
                        "names": names,
                        "title": corrected_title,
                        "journal": journal,
                        "pub_date": pub_date[0:4],
                        "volume": volume,
                        "issue": issue,
                        "pages": pages,
                        "doi": doi,
                    }
                )
                
            # create exception
            except ValueError:
                allArticles.append(
                    {
                        "pubmed_id": pubmed_id,
                        "names": "No authors found",
                        "title": title,
                        "journal": journal,
                        "pub_date": pub_date[0:4],
                        "volume": volume,
                        "issue": issue,
                        "pages": pages,
                        "doi": doi,
                    }
                )
        return jsonify({'status': 'ok', 'data': allArticles})


if __name__ == '__main__':
    app.run(host='0.0.0.0',port=5000,debug=True)
