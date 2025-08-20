
import abcEconomics as abce

class Household(abce.Agent):
    """
    가계(Household) 에이전트입니다.
    이 에이전트는 노동을 공급하고, 그 대가로 얻은 소득(임금+배당금)을 사용하여
    상품(빵)을 소비함으로써 효용을 극대화하는 것을 목표로 합니다.
    """
    def init(self, num_firms):
        """
        가계 에이전트의 초기 상태를 설정합니다.

        Args:
            num_firms (int): 시뮬레이션에 존재하는 총 기업의 수.
                             기업의 이윤을 균등하게 배당받기 위해 필요합니다.
        """
        # 모든 가계는 1 단위의 노동(labor)을 보유하며, 이를 시장에 판매할 수 있습니다.
        self.create('labor', 1)
        
        # 초기 자본으로 약간의 돈(money)을 보유합니다.
        self.create('money', 10)

        # 소비할 상품인 빵(bread)은 0개에서 시작합니다.
        self.create('bread', 0)

        # 가계는 모든 기업의 지분을 1/n 씩 균등하게 소유합니다.
        # 이 변수는 기업으로부터 배당금을 받을 때 사용됩니다.
        self.ownership_share = 1.0 / num_firms

        # 이력(log)을 기록하여 시간 경과에 따른 자산 변화를 관찰합니다.
        self.log_data = {'money': [], 'bread': [], 'utility': []}


    def sell_labor(self):
        """
        노동 시장에 나와있는 일자리 제안(job offers)을 확인하고 노동을 판매합니다.
        이 모델에서는 간단하게, 가장 높은 임금을 제시하는 기업에게 노동을 판매합니다.
        """
        # 'labor' 재화에 대한 모든 구매 오퍼(채용 공고)를 가져옵니다.
        offers = self.get_offers('labor')
        
        # 만약 오퍼가 있고, 내가 팔 수 있는 노동력이 있다면
        if offers and self.possession('labor') > 0:
            # 임금이 높은 순으로 오퍼를 정렬합니다.
            best_offer = max(offers, key=lambda offer: offer['price'])
            
            # 가장 좋은 오퍼를 수락합니다.
            self.accept(best_offer)
            # print(f"{self.name} accepted job offer from {best_offer['sender']} at wage {best_offer['price']:.2f}")


    def buy_bread(self):
        """
        상품 시장에 나와있는 빵(bread)을 구매하여 효용을 극대화합니다.
        가진 돈(money)의 한계 내에서 빵을 구매합니다.
        """
        # 'bread' 재화에 대한 모든 판매 오퍼를 가져옵니다.
        offers = self.get_offers('bread')

        # 오퍼가 있고, 내가 돈이 있다면
        if offers and self.possession('money') > 0:
            # 가격이 낮은 순으로 오퍼를 정렬합니다.
            best_offer = min(offers, key=lambda offer: offer['price'])

            # 가진 돈으로 살 수 있는 만큼 빵을 구매합니다.
            # 이 부분은 더 정교한 효용 극대화 로직으로 발전시킬 수 있습니다.
            # 예를 들어, (한계효용 / 가격)을 계산하여 구매 결정을 내릴 수 있습니다.
            # 여기서는 단순화를 위해, 가장 싼 빵을 1개 구매 시도합니다.
            quantity_to_buy = 1
            if self.possession('money') >= best_offer['price'] * quantity_to_buy:
                self.buy(best_offer, quantity_to_buy)
                # print(f"{self.name} bought bread from {best_offer['sender']} at price {best_offer['price']:.2f}")


    

    def consumption(self):
        """
        보유한 빵을 소비하여 효용을 얻고, 소비된 빵은 사라집니다.
        """
        # 이 단계는 모델의 복잡도에 따라 추가/수정될 수 있습니다.
        # 현재는 빵을 구매하면 바로 효용이 발생하는 것으로 간주하고,
        # 빵은 다음 라운드로 이월(저축)된다고 가정합니다.
        # 따라서 이 함수에서는 특별한 작업을 수행하지 않습니다.
        pass

    def log_status(self):
        """
        매 라운드 종료 시점의 상태를 기록합니다.
        """
        money = self.possession('money')
        bread = self.possession('bread')
        
        # 효용함수 U = log(1 + bread)
        # 빵이 0개일 때 log(1)=0, 빵이 증가할수록 효용은 증가하지만 증가율은 감소.
        import math
        utility = math.log(1 + bread)

        self.log_data['money'].append(money)
        self.log_data['bread'].append(bread)
        self.log_data['utility'].append(utility)
